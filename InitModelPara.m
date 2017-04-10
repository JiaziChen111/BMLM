function modelPara = InitModelPara(dimPara, hyperPara, optPara)

[N, T, J, D1, D2, NT, TJ, NTJ] = UnpackDimPara(dimPara);
[alpha, b_0beta, kappa_0beta, upsilon_0beta, S_0beta, LS_0beta, Sigma_0beta, invSigma_0beta, LSigma_0beta, ...
    b_0theta, kappa_0theta, upsilon_0theta, S_0theta, r_0theta, delta_0theta] ...
    = UnpackHyperPara(hyperPara);
[dataFile, initFile, outPath, flagMH, flagSP, flagCS, flagPar, rngSeed, n_burnin, n_collect, n_MH, prtIntv] = UnpackOptPara(optPara);

% b_theta = zeros(1, D2);
% Sigma_theta = S_0theta;

% We initialized the parameters by fitting an individual logit model.
if ~isempty(initFile)
    coef = csvread(initFile);
    beta = coef(1 : N, :);
else
    beta = randn(N, D1);
end
C = 10;
b_beta = zeros(C, D1);
[c, b_beta(:, 1 : 3)] = kmeans(beta(:, 1 : 3), C);
Sigma_beta = repmat(Sigma_0beta, [1, 1, C]);

% theta = repmat(b_theta, N, 1); 


if flagSP == 1
    % lambda ~ Gamma(r_0theta, delta_0theta) -- Note: Matlab's gamma is
    % different from that in Bayesian LASSO
    lambda_theta = random('gam', r_0theta, 1 / delta_0theta, 1);
    mu = 1 / lambda_theta;
    tau_theta = random('exp', repmat(mu, [1, D2]));
    tau_theta = min(tau_theta, 1e3);
    Sigma_theta = diag(tau_theta) / kappa_0theta;
    b_theta = zeros(1, D2);
    theta = mvnrnd(repmat(b_theta, [N, 1]), Sigma_theta);
end

if flagSP == 0
    % lambda ~ Gamma(r_0theta, delta_0theta) -- Note: Matlab's gamma is
    % different from that in Bayesian LASSO
    if flagCS == 0 % each customer has its own lambda and set of tau
        lambda_theta = random('gam', r_0theta, 1 / delta_0theta, [N, 1]);
        % tau ~ exponential(mu)
        mu = 2 ./ lambda_theta;
        %tau_theta = random('exp', mu);
        tau_theta = random('exp', repmat(mu, [1, D2]));
        tau_theta = min(tau_theta, 1e3);
        Sigma_theta = [];
        b_theta = [];
        theta = randn(N, D2) .* sqrt(tau_theta);
    end
    
%     if flagCS == 0 % each product has its own lambda and set of tau
%         lambda_theta = ones(1, D2); %random('gam', r_0theta, 1 / delta_0theta, [1, D2]);
%         % tau ~ exponential(mu)
%         mu = 1 ./ lambda_theta; %2 ./ lambda_theta;
%         %tau_theta = random('exp', mu);
%         tau_theta = random('exp', repmat(mu, [N, 1]));
%         tau_theta = min(tau_theta, 1e3);
%         Sigma_theta = [];
%         b_theta = [];
%         theta = randn(N, D2) .* sqrt(tau_theta);
%     end
    
    if flagCS == 1 % one lambda, one set of tau for all customers
        lambda_theta = 1;% random('gam', r_0theta, 1 / delta_0theta, [1, 1]);%r_0theta / delta_0theta; 
        % tau ~ exponential(mu)
        mu = 1 / lambda_theta; %2 / lambda_theta;
        %tau_theta = random('exp', mu);
        tau_theta = random('exp', repmat(mu, [1, D2]));
        %tau_theta = min(tau_theta, 1e3);
        tau_theta = repmat(tau_theta, [N, 1]); % replicate
        Sigma_theta = [];
        b_theta = [];
        theta = randn(N, D2) .* sqrt(tau_theta);
    end
end

if flagSP == -1
    b_theta = zeros(1, D2);
    Sigma_theta = S_0theta;
    theta = randn(N, D2); %repmat(b_theta, N, 1);
    tau_theta = [];
    lambda_theta = [];
end

% modelPara = PackModelPara(c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta);
modelPara = PackModelPara(c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta);

return;