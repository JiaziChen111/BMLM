function [modelPara, compPara, pcptPara] = UpdateBetaPrior(dimPara, hyperPara, optPara, modelPara, compPara, pcptPara)

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

%% Update b_theta and Sigma_theta
%fprintf('number of clusters = %d\n', C);
for k = 1 : C    
    if (1) % learn Sigma_beta
        % 1.c Update Sigma_beta
        % Note: Matlab's IW is different from that in Train (2003). Need to adjust parameters accordingly.
        dfPost = upsilon_0beta + compSize(k);
        meanBeta = mean(beta(find(c==k), :), 1);
        vec = beta(find(c==k), :) - repmat(meanBeta, compSize(k), 1);
        vec2 = meanBeta - b_0beta;
        scalePost = S_0beta + vec' * vec + ...
            kappa_0beta * compSize(k) / (kappa_0beta + compSize(k)) * (vec2' * vec2);
        % draw a new Sigma_beta for the current component
        Sigma_beta(:, :, k) = iwishrnd(scalePost, dfPost);
        % 1.b Update b_beta
        kappaPost = kappa_0beta + compSize(k);
        meanPost = (kappa_0beta * b_0beta + compSize(k) * meanBeta) / (kappa_0beta + compSize(k));
        b_beta(k, :) = mvnrnd(meanPost, Sigma_beta(:, :, k) / kappaPost);
        
    else % fix Sigma_beta, to prevent large value of beta/b_beta
        Sigma_beta(:, :, k) = eye(D1);
        meanBeta = mean(beta(find(c==k), :), 1);
        covPost = inv(inv(S_0beta) + compSize(k) * inv(Sigma_beta(:, :, k)));
        meanPost = covPost * (inv(S_0beta) * b_0beta' + compSize(k) * inv(Sigma_beta(:, :, k)) * meanBeta');
        b_beta(k, :) = mvnrnd(meanPost', covPost);
    end
    
end

%% Repack parameters
modelPara = PackModelPara(c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta);
compPara = PackCompPara(C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta);
pcptPara = PackPcptPara(detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
    detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL);


