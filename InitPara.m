clear; clc; close all;
delete(gcp);
parpool;
ms.UseParallel = true;

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
% c:                Length-N column vector
% Hyper-parameters:
% b_beta:           C X D1 matrix
% Sigma_beta:       D1 X D1 X C matrix
% alpha:            Non-negative scalar
% b_0beta:          Length-D1 row vector
% Sigma_0beta:      D1 X D1 matrix
% upsilon_0beta:    Scalar
% S_0beta:          D1 X D1 matrix
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
D1 = 2;
D2 = 0; % no store factors in this model


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

% Don't add X2 to the model
X2 = [];
% % X2
% X2 = eye(J - 1); % identity matrix of dim J
% X2 = [X2; zeros(1, J - 1)]; % "Other" stores

% Normalize X1
% the following will generate an error in  parallel computing
% parfor i = 1 : N
%     for t = 1 : T
%         X1(i, t, :, :) = X1(i, t, :, :) - repmat(X1(i, t, J, :), 1, 1, J, 1);
%     end
% end
% use this instead:
X1New = zeros(J, D1, T, N);
parfor i = 1 : N
    for t = 1 : T
        X1New(:, :,  t, i) = X1(:, :, t, i) - repmat(X1(J, :, t, i), J, 1);
    end
end
X1 = X1New;
 

% In Martin's code, the initial value of beta and theta is loaded from a
% file. We don't have it so we initialize the parameters by fitting a logit
% model.
coef = zeros(N, D1 + D2);
parfor i = 1: N
    [coef(i, :), ~] = FitLogit(X1(:, :, :, i), X2, Q(:, :, i));
end

csvwrite(sprintf('ParaInitValue_%dVars.csv', D1), coef);

delete(gcp);
