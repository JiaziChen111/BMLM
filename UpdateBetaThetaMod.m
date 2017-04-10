function [modelPara, compPara, pcptPara] = UpdateBetaThetaMod(dimPara, hyperPara, optPara, modelPara, compPara, pcptPara, X1S, X2, Q)

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

%% Set parameters for this step
% MH factors
rho_beta = 0.5;
rho_theta = 0.1;
% PG truncation level
trunc = 2000;

%% Update beta and theta
if flagMH == 1 % MH on Train (2003), Page 304
    % 2.a Update beta
    % The proposal distribution is N(beta, RHO^2*Sigma_beta)
    % Ya's implementation
    for iterMH = 1 : n_MH
        betaNew = MyMvnrnd(beta, rho_beta^2 * Sigma_beta(:, :, c), rho_beta * LSigma_beta(:, : , c));
        MNLLLNew = MNLLikelihood(X1S, prod2, Q, betaNew, [], 1, 2);
        logLLTerm = MNLLLNew + ...
            NormalLikelihoodBeta(betaNew,  b_beta(c, :), [], invSigma_beta(:, :, c), detSigma_beta(c), 2, 1, 1) - ...
            MNLLL - ...
            NLLBeta;
        % decide to accept or reject the new draw
        pAcc = exp(logLLTerm);
        flagAcc = AcceptReject(pAcc);
        idx = find(flagAcc > 0);
        beta(idx, :) = betaNew(idx, :);
        MNLLL(idx, :) = MNLLLNew(idx, :);
        if length(idx) > 1
            NLLBeta(idx) = NormalLikelihoodBeta(beta(idx, :),  b_beta(c(idx), :), [], invSigma_beta(:, :, c(idx)), detSigma_beta(c(idx)), 2, 1, 1);
        end
    end
    % 2.b Update theta
    for iterMH = 1 : n_MH
        thetaNew = mvnrnd(theta, rho_theta^2 * Sigma_theta);
        logLLTerm = MNLLikelihood(prod1, X2, Q, [], thetaNew, 1, 1) + ...
            NormalLikelihoodTheta(thetaNew,  b_theta, Sigma_theta, detSigma_theta, 1) - ...
            MNLLL - ...
            NormalLikelihoodTheta(theta,  b_theta, Sigma_theta, detSigma_theta, 1);
        % decide to accept or reject the new draw
        pAcc = exp(logLLTerm);
        flagAcc = AcceptReject(pAcc);
        theta(flagAcc > 0, :) = thetaNew(flagAcc > 0, :);
    end
end

% Polya-Gamma new MH algorithm by Ya
% The proposal distribution is a polya-gamma distribution
% A set of (beta, theta) for each (i)
% A set of (omega) for each (i, t)

% %% old fast
% if flagMH == 0 &&(0)
%     etaTrue = prod1(1 : (J - 1), :, :) + prod2(1 : (J - 1), :, :);
%     etaProp = squeeze(sum(qNorm .* permute(etaTrue, [2, 1, 3]), 2)); %weighted sum, T X N matrix
%     betaNew = zeros(N, D1);
%     thetaNew = zeros(N, D2);
%     for i = 1 : N
%         for t = 1 : T
%             omega(:, :, t, i) = PolyaGamRnd(qSum(t, i), etaProp(t, i), trunc);
%         end
%     end
%     invV_beta = invSigma_beta(:, : , c) + squeeze(sum(repmat(omega, [D1, D1, 1, 1]) .* X1TimesX1, 3));
%     if flagSP == 1
%         invV_theta = repmat(invSigma_theta, [1, 1, N]) + squeeze(sum(repmat(omega, [D2, D2, 1, 1]) .* X2TimesX2, 3));
%     else
%         for i = 1 : N
%             invSigma_theta(:, :, i) = diag(invTau_theta(i, :));
%         end
%         invV_theta = invSigma_theta + squeeze(sum(repmat(omega, [D2, D2, 1, 1]) .* X2TimesX2, 3));
%         invSigma_theta = [];
%     end
%     qDiff = qNJ - qSum / 2; % T X N
%     mBetaNomi = squeeze(sum(X1AvgQNorm .* repmat(reshape(qDiff, [T, 1, N]), [1, D1, 1]), 1))' + invSB_beta(c, :); % N X D1
%     if flagSP == 1
%         mThetaNomi = squeeze(sum(X2AvgQNorm .* repmat(reshape(qDiff, [T, 1, N]), [1, D2, 1]), 1))' + repmat(invSB_theta', [N, 1]); % N X D2
%     else
%         mThetaNomi = squeeze(sum(X2AvgQNorm .* repmat(reshape(qDiff, [T, 1, N]), [1, D2, 1]), 1))';       % prior is zero mean
%     end
%     for i = 1 : N
%         % beta
%         m_beta  = (invV_beta(:, :, i) \ mBetaNomi(i, :)')';
%         betaNew(i, :) = MyMvnrndPrec(m_beta, invV_beta(:, :, i));
%         % theta
%         m_theta = (invV_theta(:, :, i) \  mThetaNomi(i, :)')';
%         thetaNew(i, :) = MyMvnrndPrec(m_theta, invV_theta(:, :, i));
%     end
%     
%     % Compute the A/R rate
%     [prod1New, prod2New] = PrecProd(X1S, X2, betaNew, thetaNew, N, T, J);
%     etaTrueNew = prod1New(1 : (J - 1), :, :) + prod2New(1 : (J - 1), :, :);
%     etaPropNew = squeeze(sum(qNorm .* permute(etaTrueNew, [2, 1, 3]), 2)); %weighted sum, T X N matrix
%     pAcc = sum(qSum .* ( ...
%         log(1 +squeeze(sum(exp(etaTrue), 1))) - ...
%         log(1 +squeeze(sum(exp(etaTrueNew), 1))) + ...
%         log(1 + exp(etaPropNew)) - ...
%         log(1 + exp(etaProp))), 1)'; % N X 1 vector
%     pAcc = min(exp(pAcc), 1);
%     flagAcc = AcceptReject(pAcc);
%     idx = find(flagAcc > 0);
%     beta(idx, :) = betaNew(idx, :);
%     theta(idx, :) = thetaNew(idx, :);
%     % Update
%     if length(idx) > 1
%         NLLBeta(idx) = NormalLikelihoodBeta(beta(idx, :),  b_beta(c(idx), :), [], invSigma_beta(:, :, c(idx)), detSigma_beta(c(idx)), 2, 1, 1);
%     end
% end
% 
% %% old
% if flagMH == 0 && (0)
%     etaTrue = prod1(1 : (J - 1), :, :) + prod2(1 : (J - 1), :, :);
%     etaTrueNew = zeros(J - 1, T, N);
%     etaProp = zeros(T, N);
%     etaPropNew = zeros(T, N);
%     betaNew = zeros(N, D1);
%     thetaNew = zeros(N, D2);
%     for i = 1 : N
%         omega = zeros(1, T);
%         invV_beta = invSigma_beta(:, : , c(i));
%         invV_theta = invSigma_theta;
%         for t = 1 : T
%             etaProp(t, i) = qNorm(t, :, i) * etaTrue(:, t, i); %weighted sum
%             omega(1, t) = PolyaGamRnd(qSum(t, i), etaProp(t, i), trunc);
%             invV_beta  = invV_beta + omega(1, t) * X1TimesX1(:, :, t, i);
%             invV_theta  = invV_theta + omega(1, t) * X2TimesX2(:, :, t, i);
%         end
%         % beta
%         m_beta  = (invV_beta \ (X1AvgQNorm(:, :, i)' * (qNJ(:, i) - qSum(:, i) / 2) + invSB_beta(c(i), :)'))';
%         betaNew(i, :) = MyMvnrndPrec(m_beta, invV_beta);
%         % theta
%         m_theta = (invV_theta \ (X2AvgQNorm(:, :, i)' * (qNJ(:, i) - qSum(:, i) / 2) + invSB_theta))';       % prior is zero mean
%         thetaNew(i, :) = MyMvnrndPrec(m_theta, invV_theta);
%     end
%     % Compute the A/R rate
%     [prod1New, prod2New] = PrecProd(X1S, X2, betaNew, thetaNew, N, T, J);
%     etaTrueNew = prod1New(1 : (J - 1), :, :) + prod2New(1 : (J - 1), :, :);
%     pAcc = zeros(N, 1);
%     for i = 1 : N
%         for t = 1 : T
%             etaPropNew(t, i) = qNorm(t, :, i) * etaTrueNew(:, t, i); %weighted sum
%             pAcc(i, 1)  = pAcc(i, 1) + qSum(t, i) * ( ...
%                 log(1 + sum(exp(etaTrue(:, t, i)))) - ...
%                 log(1 + sum(exp(etaTrueNew(:, t, i)))) + ...
%                 log(1 + exp(etaPropNew(t, i))) - ...
%                 log(1 + exp(etaProp(t, i))));
%         end
%         pAcc(i, 1) = min(exp(pAcc(i, 1)), 1);
%     end
%     flagAcc = AcceptReject(pAcc);
%     idx = find(flagAcc > 0);
%     beta(idx, :) = betaNew(idx, :);
%     theta(idx, :) = thetaNew(idx, :);
%     % Update
%     if length(idx) > 1
%         NLLBeta(idx) = NormalLikelihoodBeta(beta(idx, :),  b_beta(c(idx), :), [], invSigma_beta(:, :, c(idx)), detSigma_beta(c(idx)), 2, 1, 1);
%     end
% end

% %% new
% if flagMH == 0
%     etaTrue = prod1(1 : (J - 1), :, :) + prod2(1 : (J - 1), :, :);
%     etaTrueNew = zeros(J - 1, T, N);
%     etaProp = zeros(T, N);
%     etaPropNew = zeros(T, N);
%     betaNew = zeros(N, D1);
%     thetaNew = zeros(N, D2);
%     % update beta
%     for i = 1 : N
%         omega = zeros(1, T);
%         invV_beta = invSigma_beta(:, : , c(i));
%         for t = 1 : T
%             etaProp(t, i) = qNorm(t, :, i) * etaTrue(:, t, i); %weighted sum
%             omega(1, t) = PolyaGamRnd(qSum(t, i), etaProp(t, i), trunc);
%             invV_beta  = invV_beta + omega(1, t) * X1TimesX1(:, :, t, i);
%         end
%         % beta
%         m_beta  = (invV_beta \ ...
%             (X1AvgQNorm(:, :, i)' * (qNJ(:, i) - qSum(:, i) / 2) + ...
%             X1AvgQNorm(:, :, i)' * diag(omega) * X2AvgQNorm(:, :, i) * theta(i, :)' + ...
%             invSB_beta(c(i), :)') ...
%             )';
%         betaNew(i, :) = MyMvnrndPrec(m_beta, invV_beta);
%     end
%     % Compute the A/R rate
%     [prod1New, prod2New] = PrecProd(X1S, X2, betaNew, theta, N, T, J);
%     etaTrueNew = prod1New(1 : (J - 1), :, :) + prod2(1 : (J - 1), :, :);
%     pAcc = zeros(N, 1);
%     for i = 1 : N
%         for t = 1 : T
%             etaPropNew(t, i) = qNorm(t, :, i) * etaTrueNew(:, t, i); %weighted sum
%             pAcc(i, 1)  = pAcc(i, 1) + qSum(t, i) * ( ...
%                 log(1 + sum(exp(etaTrue(:, t, i)))) - ...
%                 log(1 + sum(exp(etaTrueNew(:, t, i)))) + ...
%                 log(1 + exp(etaPropNew(t, i))) - ...
%                 log(1 + exp(etaProp(t, i))));
% %             log(1 + sum(exp(etaTrue(:, t, i))))
% %              log(1 + sum(exp(etaTrueNew(:, t, i)))) 
% %              log(1 + exp(etaPropNew(t, i))) 
% %                 log(1 + exp(etaProp(t, i)))
%         end
%         pAcc(i, 1) = min(exp(pAcc(i, 1)), 1);
%     end
%     flagAcc = AcceptReject(pAcc);
%     idx = find(flagAcc > 0);
% %     fprintf(['accepted ', num2str(length(idx)), '\n']);
%     beta(idx, :) = betaNew(idx, :);
%     etaTrue(:, :, idx) = etaTrueNew(:, :, idx);
%     % Update
%     if length(idx) > 1
%         NLLBeta(idx) = NormalLikelihoodBeta(beta(idx, :),  b_beta(c(idx), :), [], invSigma_beta(:, :, c(idx)), detSigma_beta(c(idx)), 2, 1, 1);
%     end
%     
%     % update theta
%      for i = 1 : N
%          if flagSP == 0
%              invV_theta = invSigma_theta(:, :, i);
%          else
%              invV_theta = invSigma_theta;
%          end
%         for t = 1 : T
%             etaProp(t, i) = qNorm(t, :, i) * etaTrue(:, t, i); %weighted sum
%             omega(1, t) = PolyaGamRnd(qSum(t, i), etaProp(t, i), trunc);
%             invV_theta  = invV_theta + omega(1, t) * X2TimesX2(:, :, t, i);
%         end
%         m_theta(i, :) = (invV_theta \ ...
%             (X2AvgQNorm(:, :, i)' * (qNJ(:, i) - qSum(:, i) / 2) + ...
%             X2AvgQNorm(:, :, i)' * diag(omega) * X1AvgQNorm(:, :, i) * betaNew(i, :)' + ...
%             invSB_theta))';       
%         thetaNew(i, :) = MyMvnrndPrec(m_theta(i, :), invV_theta);
%     end
%     % Compute the A/R rate
%     [prod1New, prod2New] = PrecProd(X1S, X2, beta, thetaNew, N, T, J);
%     etaTrueNew = prod1New(1 : (J - 1), :, :) + prod2New(1 : (J - 1), :, :);
%     pAcc = zeros(N, 1);
%     for i = 1 : N
%         for t = 1 : T
%             etaPropNew(t, i) = qNorm(t, :, i) * etaTrueNew(:, t, i); %weighted sum
%             pAcc(i, 1)  = pAcc(i, 1) + qSum(t, i) * ( ...
%                 log(1 + sum(exp(etaTrue(:, t, i)))) - ...
%                 log(1 + sum(exp(etaTrueNew(:, t, i)))) + ...
%                 log(1 + exp(etaPropNew(t, i))) - ...
%                 log(1 + exp(etaProp(t, i))));
% %             log(1 + sum(exp(etaTrue(:, t, i))))
% %              log(1 + sum(exp(etaTrueNew(:, t, i)))) 
% %              log(1 + exp(etaPropNew(t, i))) 
% %                 log(1 + exp(etaProp(t, i)))
%         end
%         pAcc(i, 1) = min(exp(pAcc(i, 1)), 1);
%     end
%     flagAcc = AcceptReject(pAcc);
%     idx = find(flagAcc > 0);
% %     fprintf(['accepted ', num2str(length(idx)), '\n']);
%     theta(idx, :) = thetaNew(idx, :);
% end


%% new fast
if flagMH == 0
    betaNew = zeros(N, D1);
    thetaNew = zeros(N, D2);
    qDiff = qNJ - qSum / 2; % T X N
    
    % Update beata
    etaTrue = prod1(1 : (J - 1), :, :) + prod2(1 : (J - 1), :, :);
    etaProp = squeeze(sum(qNorm .* permute(etaTrue, [2, 1, 3]), 2)); %weighted sum, T X N matrix
    for i = 1 : N
        for t = 1 : T
            omega(:, :, t, i) = PolyaGamRnd(qSum(t, i), etaProp(t, i), trunc);
        end
    end
    invV_beta = invSigma_beta(:, : , c) + squeeze(sum(repmat(omega, [D1, D1, 1, 1]) .* X1TimesX1, 3));    
    X1omega = X1AvgQNorm .* repmat(reshape(omega, [T, 1, N]), [1, D1, 1]);
    X2theta = sum(X2AvgQNorm .* repmat(reshape(theta, [1, J-1, N]), [T, 1, 1]),2);
    ProdX1X2 = squeeze(sum(X1omega.*repmat(X2theta,[1,D1,1]),1))';
    mBetaNomi = squeeze(sum(X1AvgQNorm .* repmat(reshape(qDiff, [T, 1, N]), [1, D1, 1]), 1))' - ProdX1X2 + invSB_beta(c, :); % N X D1
    for i = 1 : N
        m_beta  = (invV_beta(:, :, i) \ mBetaNomi(i, :)')';
        betaNew(i, :) = MyMvnrndPrec(m_beta, invV_beta(:, :, i));
    end
    % Compute the A/R rate
    [prod1New, ~] = PrecProd(X1S, X2, betaNew, theta, N, T, J);
    etaTrueNew = prod1New(1 : (J - 1), :, :) + prod2(1 : (J - 1), :, :);
    etaPropNew = squeeze(sum(qNorm .* permute(etaTrueNew, [2, 1, 3]), 2)); %weighted sum, T X N matrix
    pAcc = sum(qSum .* ( ...
        log(1 +squeeze(sum(exp(etaTrue), 1))) - ...
        log(1 +squeeze(sum(exp(etaTrueNew), 1))) + ...
        log(1 + exp(etaPropNew)) - ...
        log(1 + exp(etaProp))), 1)'; % N X 1 vector
    pAcc = min(exp(pAcc), 1);
    flagAcc = AcceptReject(pAcc);
    idx = find(flagAcc > 0);
    %     fprintf(['accepted ', num2str(length(idx)), '\n']);
    beta(idx, :) = betaNew(idx, :);
    etaTrue(:, :, idx) = etaTrueNew(:, :, idx); %update
    etaProp(:, idx) = etaPropNew(:, idx); %update
    % Update
    if length(idx) > 1
        NLLBeta(idx) = NormalLikelihoodBeta(beta(idx, :),  b_beta(c(idx), :), [], invSigma_beta(:, :, c(idx)), detSigma_beta(c(idx)), 2, 1, 1);
    end
    
    % update theta
    for i = 1 : N
        for t = 1 : T
            omega(:, :, t, i) = PolyaGamRnd(qSum(t, i), etaProp(t, i), trunc);
        end
    end   
    X2omega = X2AvgQNorm .* repmat(reshape(omega, [T, 1, N]), [1, J-1, 1]);
    X1beta = sum(X1AvgQNorm .* repmat(reshape(beta, [1, D1, N]), [T, 1, 1]),2);
    ProdX1X2 = squeeze(sum(X2omega.*repmat(X1beta,[1,J-1,1]),1))';
    if flagSP == 0
        invV_theta = invSigma_theta + squeeze(sum(repmat(omega, [D2, D2, 1, 1]) .* X2TimesX2, 3));
    else
        invV_theta = repmat(invSigma_theta, [1, 1, N]) + squeeze(sum(repmat(omega, [D2, D2, 1, 1]) .* X2TimesX2, 3));
    end
    if flagSP == 1
        mThetaNomi = squeeze(sum(X2AvgQNorm .* repmat(reshape(qDiff, [T, 1, N]), [1, D2, 1]), 1))' -ProdX1X2 + repmat(invSB_theta', [N, 1]); % N X D2
    else
        mThetaNomi = squeeze(sum(X2AvgQNorm .* repmat(reshape(qDiff, [T, 1, N]), [1, D2, 1]), 1))'-ProdX1X2;       % prior is zero mean
    end
    for i = 1 : N
        m_theta = (invV_theta(:, :, i) \  mThetaNomi(i, :)')';
        thetaNew(i, :) = MyMvnrndPrec(m_theta, invV_theta(:, :, i));
    end  
    % Compute the A/R rate
    [prod1New, prod2New] = PrecProd(X1S, X2, beta, thetaNew, N, T, J);
    etaTrueNew = prod1New(1 : (J - 1), :, :) + prod2New(1 : (J - 1), :, :);
    etaPropNew = squeeze(sum(qNorm .* permute(etaTrueNew, [2, 1, 3]), 2)); %weighted sum, T X N matrix
    pAcc = sum(qSum .* ( ...
        log(1 +squeeze(sum(exp(etaTrue), 1))) - ...
        log(1 +squeeze(sum(exp(etaTrueNew), 1))) + ...
        log(1 + exp(etaPropNew)) - ...
        log(1 + exp(etaProp))), 1)'; % N X 1 vector
    flagAcc = AcceptReject(pAcc);
    idx = find(flagAcc > 0);
    %     fprintf(['accepted ', num2str(length(idx)), '\n']);
    theta(idx, :) = thetaNew(idx, :);
end


%% Repack parameters
modelPara = PackModelPara(c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta);
compPara = PackCompPara(C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta);
pcptPara = PackPcptPara(detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
    detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL);


