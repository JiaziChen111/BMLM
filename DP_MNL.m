function DP_MNL(flagPar, flagCoda, rngSeed, flagSplitBT, flagNIW, n_burnin, n_collect, n_MH, D1, outPath)
% Dirichlet Process Mixture (DPM) + Multinomial Logit (MNL)
% All the input paramters must be specified. They could be empty though.

if nargin < 9
    fprintf('Error: All the input paramters must be specified. They could be empty though.\n');
    return;
end

if isempty(flagPar)
    flagPar = true;
end

if isempty(flagCoda)
    flagCoda = false;
end

if isempty(rngSeed)
    rngSeed = 321654987;
end

% Update beta and theta seperately
flagSplitBT = true;

% use NIW prior on theta
flagNIW = true;

if isempty(n_burnin)
    n_burnin = 5000;
end

if isempty(n_collect)
    n_collect = 5000;
end

if isempty(n_MH)
    n_MH = 50;
end

if isempty(D1)
    D1 = 2;
end

if isempty(outPath)
    outPath = sprintf('nBurnin_%d_nCollect_%d_nMH_%d_splitBT_%d_NIW_%d_D1_%d', n_burnin, n_collect, n_MH, flagSplitBT, flagNIW, D1);
    outPath = fullfile('Output', outPath);
    if exist(outPath) ~= 7
        mkdir(outPath);
    end
end

if flagPar
    delete(gcp);
    parpool;
    ms.UseParallel = true;
    %myCluster = parcluster('local'); 
    %delete(myCluster.Jobs);
end

rng(rngSeed); % It may not generate the same results even if we set the seed for the random generator the same as in Martin's code.

flagCheckError = false;

% Notations follow Matt's paper
% N:                Number of individuals
% T:                Number of time periods
% J:                Number of choices
% D1:               Number of observed variables
% D2:               Number of choice-specific variables
% C:                Number of components in DPM, i.e. clusters

% Data type:
% Observed variables:
% X_1:              J X D1 X T X N matrix (in Matlab the 1st and 2nd subscripts refers rwo and column)
% X_2:              J X D2 matrix
% Q:                T X J X N matrix
% Parameters:
% beta:             N X D1 matrix
% theta:            N X D2 matrix
% c:                Length-N column vector. beta_i = b_beta_{c_i}
% Hyper-parameters:
% b_beta:           C X D1 matrix
% alpha:            Non-negative scalar
% b_0beta:          Length-D1 row vector
% Sigma_0beta:      D1 X D1 matrix
% b_theta:          Length-D2 row vector
% Sigma_theta:      D2 X D2 matrix
% b_0theta:         Length-D2 row vector
% Sigma_0theta:     D2 X D2 matrix
% upsilon_0theta:   Scalar
% S_0theta:         D2 X D2 matrix
% Gibbs sampling parameters:
% n_burnin:         Scalar
% n_collect:        Scalar

N = 675;
T = 24; 
J = 6;
NT = N * T;
TJ = T * J;
NTJ = N * T * J;
%D1 = 2;
D2 = J - 1;


% Load data
data = load('trips_Houston_0405_pd_no_headline.txt');
if N ~=  size(data,1) / TJ
    fprintf('Error: input data size is not correct!');
end

% Q: #visits to a store
for i = 1 : N
    Q(:, :, i) = reshape(data((i - 1) * TJ + (1 : TJ), 5), J, T)';
end

% X1: 
X1 = zeros(J, D1, T, N);
for i = 1 : N
    for t = 1 : T
        X1(:, 1, t, i) = prod(data((i - 1) * TJ + (t - 1) * J + (1 : J), 7 : 8), 2); % price*distance
        if D1 > 1
            X1(:, 2, t, i) = data((i - 1) * TJ + (t - 1) * J + (1 : J), 8); %distance
        end
        if D1 > 2
            X1(:, 3, t, i) = data((i - 1) * TJ + (t - 1) * J + (1 : J), 7); %price
        end
    end
end

% X2
X2 = eye(J - 1); % identity matrix of dim J
X2 = [X2; zeros(1, J - 1)]; % "Other" stores

% Normalize X1
X1New = zeros(J, D1, T, N);
for i = 1 : N
    for t = 1 : T
        X1New(:, :,  t, i) = X1(:, :, t, i) - repmat(X1(J, :, t, i), J, 1);
    end
end
X1 = X1New;


% Set running options:
% gibbs sampling parameters
n_total = n_burnin + n_collect;
% MH factors
rho_beta = 0.5;
rho_theta = 0.1;

% Set hyper-parameters
alpha = 1e-3;
b_0beta = zeros(1.0, D1);
Sigma_0beta = eye(D1) * 25.0;
invSigma_0beta = inv(Sigma_0beta);
% for theta
if flagNIW % using normal-inverse-wishart
    kappa_0theta = 1;    
else % Martin's
    Sigma_0theta = eye(D2) * 10;
    invSigma_0theta = inv(Sigma_0theta);
end
b_0theta = zeros(1, D2);
upsilon_0theta = D2 + 1;
S_0theta = eye(D2) * upsilon_0theta; % Note: Matlab's IW is different from that in Train (2003). Need to adjust parameters accordingly.


% Initialize parameters
b_theta = [-1.6, 0.5, -2.2, -0.3, -3.5];
Sigma_theta = S_0theta;

% for the univariate case beta and theta are initialized as a randome draw
% from a starndard normal dist and a random assignment to 10 classes
if D1 == 1
    theta = mvnrnd(zeros(N, D2), eye(D2));
    C = 10;
    b_beta = randn(C, 1);
    c = ceil(rand(N, 1) * C);
end

if D1 == 2
    % In Martin's code, the initial value of beta and theta is loaded from a
    % file. We don't have the file so we initialized the parameters by fitting a logit
    % model.
    coef = csvread(sprintf('ParaInitValue_%dVars.csv', D1));
    beta = coef;
    % Specify component parameters
    C = 10;
    b_beta(1,1) = -9.0;
    b_beta(2,1) = -7.0;
    b_beta(3,1) = -5.0;
    b_beta(4,1) = -3.0;
    b_beta(5,1) = -1.0;
    b_beta(6,1) =  1.0;
    b_beta(7,1) =  3.0;
    b_beta(8,1) =  5.0;
    b_beta(9,1) =  7.0;
    b_beta(10,1) = 9.0;
    b_beta(:,2:D1) = 0.0;
    %
    % Assign members
    c = zeros(N, 1);
    for i = 1 : N
        if (beta(i,1)<-8.0) c(i) = 1; end
        if ((beta(i,1)>-8.0) && (beta(i,1)<=-6.0)) c(i) = 2; end
        if ((beta(i,1)>-6.0) && (beta(i,1)<=-4.0)) c(i) = 3; end
        if ((beta(i,1)>-4.0) && (beta(i,1)<=-2.0)) c(i) = 4; end
        if ((beta(i,1)>-2.0) && (beta(i,1)<=-0.0)) c(i) = 5; end
        if ((beta(i,1)> 0.0) && (beta(i,1)<= 2.0)) c(i) = 6; end
        if ((beta(i,1)> 2.0) && (beta(i,1)<= 4.0)) c(i) = 7; end
        if ((beta(i,1)> 4.0) && (beta(i,1)<= 6.0)) c(i) = 8; end
        if ((beta(i,1)> 6.0) && (beta(i,1)<= 8.0)) c(i) = 9; end
        if (beta(i,1)>=8.0) c(i) = 10; end
    end

end

beta = b_beta(c, :);
theta = repmat(b_theta, N, 1); % Martin's. Why do this????

% Variables to store Gibbs samples
c_collection = zeros(N, n_collect);
beta_collection = zeros(N, D1, n_collect);
theta_collection = zeros(N, D2, n_collect);
b_beta_collection = cell(n_collect, 1);
b_theta_collection = zeros(1, D2, n_collect);

% Gibbs sampling
for iter = 1 : n_total
    if mod(iter, 1) == 0
        fprintf('iteration %d of %d...\n', iter, n_total);
    end
     
    % 1.a Update c
    % find unique cs
    compSize = zeros(C, 1);
    for j = 1 : C
        compSize(j) = sum(c == j);
    end
    
    % update every c using Algorithm 7 in Neal (2000)
    % Alg. 7 Step 1 MH
    for i = 1 : N
        % if c_i is not a singleton
        if compSize(c(i)) > 1
            % draw a new cluster
            b_betaNew = mvnrnd(b_0beta, Sigma_0beta);
            betaNew = b_betaNew;
            logLLTerm = MNLLikelihood(X1(:, :, :, i), X2, Q(:, :, i), betaNew, theta(i, :), 1) - ...
                MNLLikelihood(X1(:, :, :, i), X2, Q(:, :, i), beta(i, :), theta(i, :), 1);
            % decide to accept or reject the new draw
            pAcc = min(1, exp(log(alpha) - log(N - 1) + logLLTerm));
            [~, flagAcc] = AcceptReject([], [], pAcc);
            % Add this new component if accepted
            if flagAcc == 1
                [C, compSize, b_beta, ~] = AddComponent(C, compSize, b_beta, [], b_betaNew, []);
                [C, compSize, b_beta, ~, c] = MoveMember(C, compSize, b_beta, [], i, C, c);
                beta(i, :) = betaNew;
                if flagCheckError
                    errCode = CheckMembership(C, c, compSize);
                    if errCode > 0
                        fprintf('Error: membership error %d!\n', errCode);
                        return;
                    end
                end
            end           
        % singleton   
        else 
            % draw a new c_i from c_(-i)
            c_iNew = DrawCFromOtherCs(N, compSize, c(i));
            b_betaNew = b_beta(c_iNew, :);
            betaNew = b_betaNew;
            logLLTerm = MNLLikelihood(X1(:, :, :, i), X2, Q(:, :, i), betaNew, theta(i, :), 1) - ...
                MNLLikelihood(X1(:, :, :, i), X2, Q(:, :, i), beta(i, :), theta(i, :), 1);
            % decide to accept or reject the new draw
            pAcc = min(1, exp(log(N - 1) - log(alpha) + logLLTerm));
            [~, flagAcc] = AcceptReject([], [], pAcc);           
            % if accepted, move individual i to the newly-assigned compoent and drop the component where it was before since it's a singleton 
            if flagAcc == 1
                [C, compSize, b_beta, ~, c] = MoveMember(C, compSize, b_beta, [], i, c_iNew, c);
                beta(i, :) = betaNew;
                if flagCheckError
                    errCode = CheckMembership(C, c, compSize);
                    if errCode > 0
                        fprintf('Error: membership error %d!\n', errCode);
                        return;
                    end
                end
            end            
        end
    end
    
    % Alg. 7 Step 2 Partial GS
    for i = 1 : N
        % if c_i is not a singleton
        if compSize(c(i)) > 1
            F = zeros(C, 1);
            for k = 1 : C
                b_betaNew = b_beta(k, :);
                betaNew = b_betaNew;
                F(k) = MNLLikelihood(X1(:, :, :, i), X2, Q(:, :, i), betaNew, theta(i, :), 0);
            end
            c_iNew = DrawCFromOtherCs(N, compSize, c(i), F);
            if c_iNew ~= c(i)
                [C, compSize, b_beta, ~, c] = MoveMember(C, compSize, b_beta, [], i, c_iNew, c);
                beta(i, :) = betaNew;
                if flagCheckError
                    errCode = CheckMembership(C, c, compSize);
                    if errCode > 0
                        fprintf('Error: membership error %d!\n', errCode);
                        return;
                    end
                end
            end
        end
    end
     
    % 2 Update b_beta(beta)
    % MH on Train (2003), Page 304
    % The proposal distribution is N(0, RHO^2*Sigma_beta)
    parfor k = 1 : C
        idxThisComp = find(c == k);
        for iterMH = 1 : n_MH
            b_betaNew = mvnrnd(b_beta(k, :), rho_beta^2 * Sigma_0beta);
            betaNew = repmat(b_betaNew, compSize(k), 1);
            logLLTerm = ...
                MNLLikelihood(X1(:, :, :, idxThisComp), X2, Q(:, :, idxThisComp), betaNew, theta(idxThisComp, :), 1) + ...
                NormalLikelihood(b_betaNew, b_0beta,  Sigma_0beta, 1) - ...
                MNLLikelihood(X1(:, :, :, idxThisComp), X2, Q(:, :, idxThisComp), beta(idxThisComp, :), theta(idxThisComp, :), 1) - ...
                NormalLikelihood(b_beta(k, :), b_0beta,  Sigma_0beta, 1);
            % decide to accept or reject the new draw
            pAcc = exp(logLLTerm);
            [~, flagAcc] = AcceptReject([], [], pAcc);
            if flagAcc == 1
                b_beta(k, :) =  b_betaNew;
            end
        end
    end
    beta = b_beta(c, :);
    
    % 2.b Update theta
    parfor i = 1 : N
        for iterMH = 1 : n_MH
            thetaNew = mvnrnd(theta(i,:), rho_theta^2 * Sigma_theta);
            logLLTerm = MNLLikelihood(X1(:, :, :, i), X2, Q(:, :, i), beta(i,:), thetaNew, 1) + ...
                NormalLikelihood(thetaNew,  b_theta, Sigma_theta, 1) - ...
                MNLLikelihood(X1(:, :, :, i), X2, Q(:, :, i), beta(i, :), theta(i, :), 1) - ...
                NormalLikelihood(theta(i, :),  b_theta, Sigma_theta, 1);
            % decide to accept or reject the new draw
            pAcc = exp(logLLTerm);
            [~, flagAcc] = AcceptReject([], [], pAcc);
            if flagAcc == 1
                theta(i,:) = thetaNew;
            end
        end
    end

    % 3&4 Update b_theta, sigma_theta
    if flagNIW % using NIW
        % 4 Update Sigma_theta
        % Note: Matlab's IW is different from that in Train (2003). Need to adjust parameters accordingly.
        dfPost = upsilon_0theta + N;
        meantheta = mean(theta, 1);
        vec = theta - repmat(meantheta, N, 1);
        vec2 = meantheta - b_0theta;
        scalePost = S_0theta + vec' * vec + ...
            kappa_0theta * N / (kappa_0theta + N) * (vec2' * vec2);
        % draw a new Sigma_theta for the current component
        Sigma_theta = iwishrnd(scalePost, dfPost);
        % 3 Update b_theta
        kappaPost = kappa_0theta + N;
        meanPost = (kappa_0theta * b_0theta + N * meantheta) / (kappa_0theta + N);
        b_theta = mvnrnd(meanPost, Sigma_theta / kappaPost);
        %
    else %Martins
        % 3 Update b_theta
        % diffuse prior on b_theta
        covPost = Sigma_theta / N;
        meanPost = mean(theta, 1);
        % draw a new b_beta for the current component
        b_theta = mvnrnd(meanPost, covPost);
        %
        % 4 Update Sigma_theta
        dfPost = upsilon_0theta + N;
        vec = theta - repmat(b_theta, N, 1);
        scalePost = S_0theta + vec'*vec;
        % draw a new Sigma_beta for the current component
        Sigma_theta = iwishrnd(scalePost, dfPost);
    end
    
    % 5 Collect samples
    if iter > n_burnin
        c_collection(:, iter - n_burnin) = c;
        beta_collection(:,:, iter - n_burnin) = beta;
        theta_collection(:,:, iter - n_burnin) = theta;
        b_beta_collection{iter - n_burnin} = b_beta; %cell
        b_theta_collection(:,:, iter - n_burnin) = b_theta;
    end
end
    
% Save samples
fileName = fullfile(outPath, 'Samples_c.mat');
save(fileName, 'c_collection');
fileName = fullfile(outPath, 'Samples_beta.mat');
save(fileName, 'beta_collection');
fileName = fullfile(outPath, 'Samples_theta.mat');
save(fileName, 'theta_collection');
fileName = fullfile(outPath, 'Samples_b_beta.mat');
save(fileName, 'b_beta_collection');
fileName = fullfile(outPath, 'Samples_b_theta.mat');
save(fileName, 'b_theta_collection');

% Evaluate convergence %need debug for univariate!!!!!!!!!!!!
if flagCoda
    addpath(genpath('Z:\Users\yx10\Toolbox\SpatialEconometrics\jplv7'));
    switch D1
        case 1
            vnames = strvcat('beta_1');
        case 2
            vnames = strvcat('beta_1','beta_2');
        case 3
            vnames = strvcat('beta_1','beta_2','beta_3');
        otherwise
            return;
    end
    info.q = 0.025;
    info.r = 0.01;
    info.s = 0.950;
    info.p1 = 0.2;
    info.p2 = 0.5;
    fid = fopen(fullfile(outPath, 'Coda_beta.txt'),'w');
    if D1 == 1
        coda(squeeze(beta_collection(1,:,:)),vnames, info, fid);
    else
        coda(squeeze(beta_collection(1,:,:))',vnames, info, fid);
    end
    fclose(fid);
    vnames = strvcat('theta_1','theta_2','theta_3','theta_4','theta_5');
    fid = fopen(fullfile(outPath, 'Coda_theta.txt'),'w');
    coda(squeeze(theta_collection(1,:,:))',vnames, info, fid);
    fclose(fid);
end


% Housekeeping
if flagPar
    delete(gcp);
end







