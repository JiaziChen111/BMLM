function initFile = FitIndividualFixLogit(fileName)

delete(gcp);
parpool(4);
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

D1 = [];
T = [];
J = [];
X1S = [];

% Load data
%fileName = 'PizzaData_N1000.mat';
load(fileName);

% Set dimensionality parameters
NT = N * T;
TJ = T * J;
NTJ = N * T * J;

% Don't add X2 to the model
X2 = [];
D2 = 0;

% In Martin's code, the initial value of beta and theta is loaded from a
% file. We don't have it so we initialize the parameters by fitting a logit
% model.
coef = zeros(N, D1 + D2);
if exist('X1', 'var')
    parfor i = 1: N
        [coef(i, :), ~] = FitLogit(X1(:, :, :, i), X2, Q(:, :, i));
    end
else
    parfor i = 1: N
        idr = (i - 1) * TJ + 1 : i * TJ;
        idc = (i - 1) * D1 + 1 : i * D1;
        xBlock = X1S(idr, idc);
        X1 = zeros(J, D1, T);
        for t = 1 : T
            X1(:, :, t) = xBlock((t-1) * J + 1 : t * J, :);
        end
        [coef(i, :), ~] = FitLogit(X1(:, :, :), X2, Q(:, :, i));
    end
end

initFile = ['ParaInitValue_', fileName(1 : end-4), '.csv'];
csvwrite(initFile, coef);

delete(gcp);
