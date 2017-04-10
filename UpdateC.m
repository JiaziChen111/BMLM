function [modelPara, compPara, pcptPara] = UpdateC(dimPara, hyperPara, optPara, modelPara, compPara, pcptPara, flagCheckError)

if ~exist('flagCheckError', 'var')
    flagCheckError = 0;
end

%% Unpack parameters
[N, T, J, D1, D2, NT, TJ, NTJ] = UnpackDimPara(dimPara);
[alpha, b_0beta, kappa_0beta, upsilon_0beta, S_0beta, LS_0beta, Sigma_0beta, invSigma_0beta, LSigma_0beta, ...
    b_0theta, kappa_0theta, upsilon_0theta, S_0theta, r_0theta, delta_0theta] ...
    = UnpackHyperPara(hyperPara);
[dataFile, initFile, outPath, flagMH, flagSP, flagCS, flagPar, rngSeed, n_burnin, n_collect, n_MH, prtIntv] = UnpackOptPara(optPara);
[c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta] = UnpackModelPara(modelPara);
if ~isempty(compPara)
    [C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta] = UnpackCompPara(compPara);
end
[detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
    detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL] = UnpackPcptPara(pcptPara);

%% Update c
if (1)
    
    % find unique cs
    compSize = zeros(C, 1);
    for j = 1 : C
        compSize(j) = sum(c == j);
    end
    compPara = PackCompPara(C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta);
    
    % update every c using Algorithm 7 in Neal (2000)
    % Alg. 7 Step 1 MH
    % First draw N new components. It's faster than drawing components one
    % by one in the loop.
    [b_betaCand, Sigma_betaCand] = DrawNewComponent(b_0beta, upsilon_0beta, kappa_0beta, S_0beta, LS_0beta, Sigma_0beta, N);
    for i = 1 : N
        % if c_i is not a singleton
        if compSize(c(i)) > 1
            % draw a new cluster
            b_betaNew = b_betaCand(i, :);
            Sigma_betaNew = Sigma_betaCand(:, :, i);
            logLLTerm = NormalLikelihoodBeta(beta(i, :), b_betaNew, Sigma_betaNew, [], [], 3, 0, 1) - ...
                NLLBeta(i); %NormalLikelihoodBeta(beta(i, :), b_beta(c(i), :), [], invSigma_beta(:, :, c(i)), detSigma_beta(c(i)), 3, 1, 1); % all the other terms are cancelled out
            % decide to accept or reject the new draw
            pAcc = min(1, alpha / (N - 1) * exp(logLLTerm));
            flagAcc = AcceptReject(pAcc);
            % Add this new component if accepted
            if flagAcc == 1
                compPara = AddComponent(compPara, b_betaNew, Sigma_betaNew);
                [C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta] = UnpackCompPara(compPara);
                [compPara, c] = MoveMember(compPara, i, C, c);
                [C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta] = UnpackCompPara(compPara);
                if flagCheckError
                    errCode = CheckMembership(C, c, compSize);
                    if errCode > 0 || sum(compSize==0) > 0
                        fprintf('Error: membership error %d!\n', errCode);
                        return;
                    end
                end
            end
            %
            % singleton
        else
            % draw a new c_i from c_(-i)
            c_iNew = DrawCFromOtherCs(N, compSize, c(i));
            b_betaNew = b_beta(c_iNew, :);
            Sigma_betaNew = Sigma_beta(:, :, c_iNew);
            logLLTerm = NormalLikelihoodBeta(beta(i, :), b_betaNew, Sigma_betaNew, [], [], 3, 0, 1) - ...
                NLLBeta(i); %NormalLikelihood(beta(i, :), b_beta(c(i), :), [], invSigma_beta(:, :, c(i)), detSigma_beta(c(i)), 3, 1, 1); % all the other terms are cancelled out
            % decide to accept or reject the new draw
            pAcc = min(1, (N - 1) / alpha * exp(logLLTerm));
            flagAcc = AcceptReject(pAcc);
            % if accepted, move individual i to the newly-assigned compoent and drop the component where it was before since it's a singleton
            if flagAcc == 1
                [compPara, c] = MoveMember(compPara, i, c_iNew, c);
                [C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta] = UnpackCompPara(compPara);
                if flagCheckError
                    errCode = CheckMembership(C, c, compSize);
                    if errCode > 0 || sum(compSize==0) > 0
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
            F = NormalLikelihoodBeta(beta(i, :), b_beta, [], invSigma_beta, detSigma_beta, 1, 1, 0);
            c_iNew = DrawCFromOtherCs(N, compSize, c(i), F);
            if c_iNew ~= c(i)
                [compPara, c] = MoveMember(compPara, i, c_iNew, c);
                [C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta] = UnpackCompPara(compPara);
                if flagCheckError
                    errCode = CheckMembership(C, c, compSize);
                    if errCode > 0 || sum(compSize==0) > 0
                        fprintf('Error: membership error %d!\n', errCode);
                        return;
                    end
                end
            end
        end
    end
end

% for debug
% c = [1, 2, 3, 4, 3, 5, 2, 3, 6, 3]';
% C = 6;
% b_beta = b_beta(1 : 6, :);
% Sigma_beta = Sigma_beta(:, :, 1 : 6);
% compSize = [1, 2, 4, 1, 1, 1];

%% Repack parameters
modelPara = PackModelPara(c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta);
compPara = PackCompPara(C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta);
pcptPara = PackPcptPara(detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
    detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL);