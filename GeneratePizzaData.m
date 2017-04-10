function fileName = GeneratePizzaData(N, J, flagPlot, sparsity)

if nargin < 3
    flagPlot = 1;
end

if nargin < 4
    sparisty = 1;
end

if J > 1e3
    delete(gcp);
    parpool(12);
    ms.UseParallel = true;
end

% Generate synthetic data to simulate frozen-pizza purchases at a store

% Assume a store sells 1000 different kinds of frozen pizza. N customers
% come buy M ([1, 10]) pizzas every day in a month (30 days). 

% A pizza is characaterized by its price, type(meat or veggi), display
% position (1 - 1000, from left to right), brand(20 brands in total).

% There are 18 latent classes in the customers based on their sensitivity
% to price (positive, negative or nuetral), preference on type (prefer
% meat, prefer veggis or nuetral), and walking distance to the
% pizza(negative or neutral). Here we assume the customers are not senstive
% to the brand.

% J:    Number of pizza choices
% N:    Number of customers
% T:    Number of periods
% D1:               Number of observed variables
% D2:               Number of choice-specific variables
% X_1:              J X D1 X T X N matrix (in Matlab the 1st and 2nd subscripts refers rwo and column)
% X_2:              J X D2 matrix
% Q:                T X J X N matrix
% C:                Number of components in DPM, i.e. latent classes
% b_beta:           C X D1 matrix
% Sigma_beta:       D1 X D1 X C matrix

filePrefix = ['PizzaData_N', num2str(N), '_J', num2str(J)];

%% Generate pizza attributes
fprintf('Generate pizza attributes...\n');
%J = 1000;
% Price: uniformly distributed between $1 and $5
pizza(:, 1) = rand(J, 1) * 4 + 1;
% Size: integer variable. 1, 2, ..., 100.
pizza(:, 2) = ceil(rand(J, 1) * 100) * 0.01; % * 0.1 to make it in the same scale with other variables
% display position: random permutation
%pizza(:, 3) = randperm(J)';
% Package quality: integer variable. 1, 2, ..., 5. USA, Italy, France, Mexico,
% Greece.
pizza(:, 4) = ceil(rand(J, 1) * 100) * 0.01;

D1 = size(pizza, 2);
D2 = J - 1;

%% Generate customer profiles
fprintf('Generate customer profiles...\n');
C = 18;
% mean coefficients
b_beta = zeros(C, D1); 
% mean coef for price
b_beta(:, 1) = [ones(1, 6) * (-4), zeros(1, 6), ones(1, 6) * 4]';
% for type
b_beta(:, 2) = repmat([-5, -5, 0, 0, 5, 5],[1, 3]); 
% for walking distance to the pizza
b_beta(:, 3) = repmat([-6, 0], [1, 9])';
% for package quality
b_beta(:, 4) = ones(C, 1) * 1;

% covariance
sigma = ones(C, D1) * 0.2; 
sigma(:, 1) = [ones(1, 6) * 0.5, ones(1, 6) * 1.5, ones(1, 6) * 0.5];
sigma(:, 2) = repmat([1, 1, 2, 2, 1, 1],[1, 3]);
sigma(:, 3) = repmat([0.5, 1], [1, 9])';
Sigma_beta = zeros(D1, D1, C);
for k = 1 : C
    Sigma_beta(:, :, k) = diag(sigma(k, :));
end

% weights of each latent class
% marginal weights
mw = {};
mw{1} = [0.3, 0.6, 0.1];
mw{2} = [0.4, 0.2, 0.4];
mw{3} = [0.2, 0.8];
% weights
w = zeros(C, 1);
iter = 1;
for i = 1 : 3
    for j = 1: 3
        for k = 1: 2
            w(iter) = mw{1}(i) * mw{2}(j) * mw{3}(k);
            iter = iter + 1;
        end
    end
end
%sum(w)

% draw beta
c = zeros(N, 1);
compSize = zeros(C, 1);
for i = 1 : N
     k = find(mnrnd(1, w) > 0); 
     c(i) = k;
     compSize(k) = compSize(k) + 1;
     beta(i, :) = mvnrnd(b_beta(k, :), Sigma_beta(:, :, k));
end

if flagPlot
    for j = 1 : 4
        [f, xi] = ksdensity(beta(:, j), 'function', 'pdf', 'width', 0.5);
        figure,
        plot(xi, f);
        title(['kde for beta_', num2str(j)]);
        xlabel(['beta_', num2str(j)]);
        ylabel('density');
        savefig([filePrefix, '_kde_beta', num2str(j), '.fig']);
        saveas(gcf, [filePrefix, '_kde_beta', num2str(j), '.jpg']);
    end
end
            
% theta
if (0)
    b_theta = zeros(1, D2);
    flag = 1;
    while flag == 1
        Sigma_theta = rand(D2, D2);
        if J <= 1e3
            Sigma_theta = (Sigma_theta + Sigma_theta') / 50;
        else
            Sigma_tsheta = (Sigma_theta + Sigma_theta') / 100;
        end
        Sigma_theta = (Sigma_theta + diag(ones(1, D2) * 1)) * 0.01;
        ei = eig(Sigma_theta);
        % check if Sigma_theta is positive definite
        if sum(ei <= 0) == 0
            flag = 0;
        end
    end
    theta = mvnrnd(b_theta, Sigma_theta, N);
else
    b_theta = zeros(1, D2);
    Sigma_theta = ones(D2, D2) * 1;
    Sigma_theta = Sigma_theta + eye(D2, D2) * 4; %variance 1
    if (0) %sparisity on product
        if sparsity < 1
            nonsparseMean = [2, 1];
            nNSM = length(nonsparseMean);
            blockSize = round(J * (1 - sparsity) / nNSM);
            for j = 1 : nNSM
                b_theta(round(J * sparsity) - 1 + ((j - 1) * blockSize + 1 : j * blockSize)) = ...
                    nonsparseMean(j);
            end
        end
        theta = mvnrnd(b_theta, Sigma_theta, N);
        % make it sparse
        theta(:, 1 : (round(J * sparsity) - 1)) = 0;
    else %sparsity on customer
        theta = mvnrnd(b_theta, Sigma_theta, N);
        % make it sparse
        theta(1 : (round(N * sparsity) - 1), :) = 0;
    end
end

%% Generate customer purchase history
fprintf('Generate customer purchase history...\n');
T = 30; % 30 days
% aisle length e.g. normalized to 0.999 so the interval between two
% adjacent pizza is 0.001
intv = 1 / J;
al = 1 - intv;
% Every time period may have a different pizza display
for t = 1 : T
    pizzaDisplay(:, t) = randperm(J)';
end

% X2
X2 = [eye(J - 1); zeros(1, J - 1)];

TJ = T * J;
X1S = sparse(N*TJ, N*D1);
Q = zeros(T, J, N);

for i = 1 : N
    if mod(i, 100) == 0
        fprintf(['created ', num2str(i), ' of ', num2str(N), ' customers'' purchase histry\n']);
    end
    
    X1 = zeros(J, D1, T);
    
    for t = 1 : T
        % price 
        X1(:, 1, t) = pizza(:, 1) + randn(J, 1) * 0.2; % price may change a little bit
        %type
        X1(:, 2, t) = pizza(:, 2);
        % walking distance to a pizza
        % the customer may start from the left or right end of the aisle.
        % Randomly chosen.
        if rand(1) < 0.5 % from the left end
            X1(:, 3, t) = (pizzaDisplay(:, t) - 1) * intv;
        else % right
            X1(:, 3, t) = al - (pizzaDisplay(:, t) - 1) * intv;
        end
        % Made in
        X1(:, 4 : end, t) = pizza(:, 4 : end);
        
        % Normalize X1 ie substract the last choice
        X1(:, :, t) = X1(:, :, t) - repmat(X1(J, :, t), J, 1);
    end
    
    % reshape X1 to improve efficiency
    xBlock = [];
    for t = 1 : T
        xBlock = [xBlock; X1(:, :, t)];
    end
    idr = (i - 1) * TJ + 1 : i * TJ;
    idc = (i - 1) * D1 + 1 : i * D1;
    X1S(idr, idc) = xBlock;
    
    % purchases
    for t = 1 : T
        % first decide the quantity (1 - 10)
        M = ceil(rand(1) * 10);       
        % then pick up the pizzas
        p = exp(X1(:, :, t) * beta(i, :)' + X2 * theta(i, :)');
        p = p / sum(p);
        pick = mnrnd(M, p);
        Q(t, :, i) = pick;
    end
end

% Choice set for each individual???
Z = zeros(N, J);
for i = 1 : N
    Z(i, :) = sum(Q(:, :, i), 1) > 0; 
end

% Pizza coverage
cvrg = sum(sum(Z, 1) > 0) / J;
fprintf('pizza coverage: %f\n', cvrg);

% Save data
fprintf('Save data...\n');
X1 = []; % save space for large data
save([filePrefix, '.mat'], 'N', 'T', 'J', 'D1', 'D2', 'X1S', 'X2', 'Q', 'Z', '-v7.3'); % only variables needed for DPM
if J < 1e3
    save([filePrefix, '_all.mat'], '-v7.3');
else
    save([filePrefix, '_all.mat'], 'pizza', 'beta', 'c', 'b_beta', 'Sigma_beta', 'theta', 'b_theta', 'Sigma_theta', '-v7.3');
end

fileName = [filePrefix, '.mat'];

delete(gcp);

return;
