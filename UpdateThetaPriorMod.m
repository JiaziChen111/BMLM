function [modelPara, compPara, pcptPara] = UpdateThetaPriorMod(dimPara, hyperPara, optPara, modelPara, compPara, pcptPara)

%% Unpack parameters
[N, T, J, D1, D2, NT, TJ, NTJ] = UnpackDimPara(dimPara);
[alpha, b_0beta, kappa_0beta, upsilon_0beta, S_0beta, LS_0beta, Sigma_0beta, invSigma_0beta, LSigma_0beta, ...
    b_0theta, kappa_0theta, upsilon_0theta, S_0theta, r_0theta, delta_0theta] ...
    = UnpackHyperPara(hyperPara);
[dataFile, initFile, outPath, flagMH, flagSP, flagCS, flagPar, rngSeed, n_burnin, n_collect, n_MH, prtIntv] = UnpackOptPara(optPara);
[c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta] = UnpackModelPara(modelPara);
[C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta] = UnpackCompPara(compPara);
[detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
    detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL] = UnpackPcptPara(pcptPara);

%%
if flagSP == 1 % sparse prior on b_theta
    % 3 Update b_theta, Sigma_theta
    dfPost = upsilon_0theta + N;
    meanTheta = mean(theta, 1);
    vec = theta - repmat(meanTheta, N, 1);
    vec2 = meanTheta; % b_0theta = all zero vector
    scalePost = diag(tau_theta) + vec' * vec + ...
        kappa_0theta * N / (kappa_0theta + N) * (vec2' * vec2);
    % draw a new Sigma_theta
    Sigma_theta = iwishrnd(scalePost, dfPost);
    kappaPost = kappa_0theta + N;
    meanPost = N * meanTheta / (kappa_0theta + N);
    b_theta = mvnrnd(meanPost, Sigma_theta / kappaPost);
    
    % Pre-compute some variables to improve efficiency
    detSigma_theta = det(Sigma_theta);
    invSigma_theta = inv(Sigma_theta);
    
    % 4 Update tau_theta (variance of theta), lambda_theta
    % Update lambda_theta ~ Gamma(r, delta)
    %     rPost = 1 + r_0theta;
    %     deltaPost = tau_theta / 2 + delta_0theta;
    rPost = D2 + r_0theta;
    deltaPost = sum(tau_theta) / 2 + delta_0theta;
    lambda_theta = random('gam', rPost, 1 / deltaPost);
    % Update tau_theta ~ inverse-gaussian(mu, lambda)
    muPost = repmat(sqrt(lambda_theta), [1, D2]) ./ abs(b_theta);
    invTau_theta = random('inversegaussian', muPost, repmat(lambda_theta, [1, D2]));
    tau_theta = 1 ./ (invTau_theta + 1e-6); % 1 X D2
end

%% sparse prior on theta
if flagSP == 0 
    % 3&4 Update tau_theta (variance of theta), lambda_theta
    if flagCS == 0 % each product has its own lambda and set of tau
        if (0)
            % each product has its own lambda and set of tau
            rPost = N + r_0theta;
            deltaPost = sum(tau_theta, 1) / 2 + delta_0theta; % D2 X 1 vector
            lambda_theta = random('gam', rPost, 1 ./ deltaPost);% D2 X 1 vector
            % Update tau_theta ~ inverse-gaussian(mu, lambda)
            muPost = repmat(sqrt(lambda_theta), [N, 1]) ./ abs(theta);
            invTau_theta = random('inversegaussian', muPost, repmat(lambda_theta, [N, 1]));
            tau_theta = 1 ./ (invTau_theta + 1e-6); % N X D2
        else
            % each customer has its own lambda and set of tau
            rPost = D2 + r_0theta;
            deltaPost = sum(tau_theta, 2) / 2 + delta_0theta; % N X 1 vector
            lambda_theta = random('gam', rPost, 1 ./ deltaPost);% N X 1 vector
            % Update tau_theta ~ inverse-gaussian(mu, lambda)
            muPost = repmat(sqrt(lambda_theta), [1, D2]) ./ abs(theta);
            invTau_theta = random('inversegaussian', muPost, repmat(lambda_theta, [1, D2]));
            tau_theta = 1 ./ (invTau_theta + 1e-6); % N X D2
        end
    end
    
    if flagCS == 1 % one lambda, one set of tau for all customers
        rPost = D2 + r_0theta;
        deltaPost = sum(mean(tau_theta, 1)) / 2 + delta_0theta; % scalar
        lambda_theta = random('gam', rPost, 1 ./ deltaPost);% 1scalar
        % Update tau_theta ~ inverse-gaussian(mu, lambda)
        muPost = repmat(sqrt(lambda_theta), [1, D2]) ./ sqrt(sum(theta.^2, 1));
        invTau_theta = random('inversegaussian', muPost, repmat(lambda_theta, [1, D2]));
        tau_theta = 1 ./ (invTau_theta + 1e-6); % 1 X D2
        tau_theta = repmat(tau_theta, [N, 1]); % replicate
    end
end

%% no sparse prior
if flagSP == -1 
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
end

%% Repack parameters
modelPara = PackModelPara(c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta);
compPara = PackCompPara(C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta);
pcptPara = PackPcptPara(detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
    detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL);


