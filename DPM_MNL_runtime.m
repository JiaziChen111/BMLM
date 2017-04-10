function [outPath,t] = DPM_MNL_runtime(optPara)
% Dirichlet Process Mixture (DPM) + Multinomial Logit (MNL)
% All the input paramters must be specified. Arguments 3 - 9 could be empty though.
%
% Notations follow Matt's paper
% N:                Number of individuals
% T:                Number of time periods
% J:                Number of choices
% D1:               Number of observed variables
% D2:               Number of choice-specific variables
% C:                Number of components in DPM, i.e. clusters
%
% Data type:
% Observed variables:
% X_1:              J X D1 X T X N matrix (in Matlab the 1st and 2nd subscripts refers rwo and column)
% X_1S:             (N*T*J) X D1 SPARSE matrix. Save space for large data. 
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
% kappa_0beta:      Scalar
% S_0beta:          D1 X D1 matrix
% upsilon_0theta:   Scalar
% kappa_0theta:     Scalar
% tau_theta:        Length-D2 row vector
% lambda_theta:     Scalar
% r_0theta:         Scalar
% delta_0theta:     Scalar
% Gibbs sampling parameters:
% n_burnin:         Scalar
% n_collect:        Scalar

%% Parsing input arguments
[dataFile, initFile, outPath, flagMH, flagSP, flagCS, flagPar, ...
    rngSeed, n_burnin, n_collect, n_MH, prtIntv] = UnpackOptPara(optPara);

if isempty(flagPar)
    flagPar = true;
end

if isempty(rngSeed)
    rngSeed = clock;
    rngSeed = round(prod(rngSeed(4:end))); %321654987;
end

if isempty(n_burnin)
    n_burnin = 5000;
end

if isempty(n_collect)
    n_collect = 5000;
end

if isempty(n_MH)
    n_MH = 1;
end

outPath = GenOutpath(optPara);

if ~exist(outPath, 'dir')
    mkdir(outPath);
end


if isempty(prtIntv)
    prtIntv = 100;
end

%% Setting up parreling computing
if flagPar
    delete(gcp);
    parpool;
    ms.UseParallel = true;
    myCluster = parcluster('local');
    delete(myCluster.Jobs);
end

%% Reset random number generator
rng(rngSeed); % It may not generate the same results even if we set the seed for the random generator the same as in Martin's code.

%% Load data
if flagPar
    N = [];
    T = [];
    J = [];
    X1 = [];
    X1S = [];
    X2 = [];
    Q = [];
    Z = [];
end
load(dataFile);
        

%% Set dimensionality parameters
NT = N * T;
TJ = T * J;
NTJ = N * T * J;
% Pack dimensionality parameters
dimPara = PackDimPara(N, T, J, D1, D2, NT, TJ, NTJ);

%% Set running options:
flagCheckError = 0;
% gibbs sampling parameters
n_total = n_burnin + n_collect;

%% Set hyper-parameters
hyperPara = SetHyperPara(dimPara, flagSP);
% [alpha, b_0beta, kappa_0beta, upsilon_0beta, S_0beta, LS_0beta, Sigma_0beta, invSigma_0beta, LSigma_0beta, ...
%     kappa_0theta, upsilon_0theta, r_0theta, delta_0theta] ...
%     = UnpackHyperPara(hyperPara);

%% Initialize model parameters
modelPara = InitModelPara(dimPara, hyperPara, optPara);
% [c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta] = UnpackModelPara(modelPara);

%% Initialize latent class parameters
compPara = [];

%% Pre-compute some variables to improve efficiency
pcptPara = Precompute(X1S, X2, Q, dimPara, modelPara, optPara, [], 1);
% [detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
%     detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL] = UnpackPcptPara(pcptPara);

%% Declare variables to store Gibbs samples
c_collection = zeros(N, n_collect);
beta_collection = zeros(N, D1, n_collect);
b_beta_collection = cell(n_collect, 1);
Sigma_beta_collection = cell(n_collect, 1);
%theta_collection = zeros(N, D2, n_collect);
% b_theta_collection = zeros(1, D2, n_collect);
% Sigma_theta_collection = cell(n_collect, 1);

%% Gibbs sampling
for iter = 1 : n_total
    if mod(iter, prtIntv) == 0
        fprintf('J = %d, N = %d, iteration %d of %d...\n', J, N, iter, n_total);
    end
    
    tic
    % 1.a Update c
    [modelPara, compPara, pcptPara] = UpdateC(dimPara, hyperPara, optPara, modelPara, compPara, pcptPara);
    t.p1 = toc;
    
    tic
    % 2. Update b_beta and Sigma_beta
    [modelPara, compPara, pcptPara] = UpdateBetaPrior(dimPara, hyperPara, optPara, modelPara, compPara, pcptPara);
    t.p2 = toc;
    
    tic
    % Pre-compute some variables to improve efficiency
    pcptPara = Precompute(X1S, X2, Q, dimPara, modelPara, optPara, pcptPara, 0);
    compPara = PackCompPara(compPara.C, compPara.compSize, compPara.b_beta, compPara.Sigma_beta, pcptPara.invSigma_beta, pcptPara.detSigma_beta);
    t.p3 = toc;
    
    tic
    % 3. Update beta and theta
    [modelPara, compPara, pcptPara] = UpdateBetaTheta(dimPara, hyperPara, optPara, modelPara, compPara, pcptPara, X1S, X2, Q);
    t.p4 = toc;
    
    tic
    % 4. Update prior of theta
    [modelPara, compPara, pcptPara] = UpdateThetaPrior(dimPara, hyperPara, optPara, modelPara, compPara, pcptPara);
    t.p5 = toc;
    
    tic
    % 5 Collect samples
    if iter > n_burnin
        c_collection(:, iter - n_burnin) = modelPara.c;
        beta_collection(:,:, iter - n_burnin) = modelPara.beta;
        b_beta_collection{iter - n_burnin} = modelPara.b_beta; %cell
        Sigma_beta_collection{iter - n_burnin} = modelPara.Sigma_beta; %cell
        if n_collect <= 100
            theta_collection(:,:, iter - n_burnin) = modelPara.theta;
            b_theta_collection(:,:, iter - n_burnin) = modelPara.b_theta;
            Sigma_theta_collection(:,:, iter - n_burnin) = modelPara.Sigma_theta;
            tau_theta_collection(:,:, iter - n_burnin) = modelPara.tau_theta;
            lambda_theta_collection(:,:, iter - n_burnin) = modelPara.lambda_theta;
        else % save space
            if mod(iter - n_burnin, ceil((iter - n_burnin) / 100)) == 1
                theta_collection(:,:, ceil((iter - n_burnin) / 100)) = modelPara.theta;
                b_theta_collection(:,:, ceil((iter - n_burnin) / 100)) = modelPara.b_theta;
                Sigma_theta_collection(:,:, ceil((iter - n_burnin) / 100)) = modelPara.Sigma_theta;
                tau_theta_collection(:,:, ceil((iter - n_burnin) / 100)) = modelPara.tau_theta;
                lambda_theta_collection(:,:, ceil((iter - n_burnin) / 100)) = modelPara.lambda_theta;
            end
        end
        
        if (mod(iter - n_burnin, 1e3) == 0 && iter - n_burnin > 0) || iter == n_burnin + n_collect 
            % Save samples
            % c
            fileName = fullfile(outPath, 'Samples_c.mat');
            save(fileName, 'c_collection');
            % beta
            fileName = fullfile(outPath, 'Samples_beta.mat');
            save(fileName, 'beta_collection');
            fileName = fullfile(outPath, 'Samples_b_beta.mat');
            save(fileName, 'b_beta_collection');
            fileName = fullfile(outPath, 'Samples_Sigma_beta.mat');
            save(fileName, 'Sigma_beta_collection');
            % theta
            fileName = fullfile(outPath, 'Samples_theta.mat');
            save(fileName, 'theta_collection');
            fileName = fullfile(outPath, 'Samples_b_theta.mat');
            save(fileName, 'b_theta_collection');
            fileName = fullfile(outPath, 'Samples_Sigma_theta.mat');
            save(fileName, 'Sigma_theta_collection');
            fileName = fullfile(outPath, 'Samples_tau_theta.mat');
            save(fileName, 'tau_theta_collection');
            fileName = fullfile(outPath, 'Samples_lambda_theta.mat');
            save(fileName, 'lambda_theta_collection');
        end
        t.p6 = toc;
    end
end


%% Housekeeping
if flagPar
    delete(gcp);
end







