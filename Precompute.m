function pcptPara = Precompute(X1S, X2, Q, dimPara, modelPara, optPara, pcptPara, flagInit)
% Pre-compute some variables to improve efficiency

[N, T, J, D1, D2, NT, TJ, NTJ] = UnpackDimPara(dimPara);
[c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta] = UnpackModelPara(modelPara);
[dataFile, initFile, outPath, flagMH, flagSP, flagCS, flagPar, rngSeed, n_burnin, n_collect, n_MH, prtIntv] = UnpackOptPara(optPara);

detSigma_beta = [];
invSigma_beta = [];
for k = 1 : C
    detSigma_beta(k, 1) = det(Sigma_beta(:, :, k));
    invSigma_beta(:, : , k) = inv(Sigma_beta(:, :, k));
    LSigma_beta(:, : , k) = chol(Sigma_beta(:, :, k), 'lower');
    invSB_beta(k, :) = b_beta(k, :) * (invSigma_beta(:, : , k))';
end
NLLBeta = NormalLikelihoodBeta(beta,  b_beta(c, :), [], invSigma_beta(:, :, c), detSigma_beta(c), 2, 1, 1);

invTau_theta = 1 ./ tau_theta;

if flagSP == 1
    invSigma_theta = inv(Sigma_theta);
    invSB_theta = Sigma_theta \ b_theta'; % D2 X 1 vector
end

if flagSP == 0
    for i = 1 : N
        invSigma_theta(:, :, i) = diag(invTau_theta(i, :));
    end
    invSB_theta = zeros(D2, 1);
end

if flagSP == -1
    invSigma_theta = inv(Sigma_theta);
    invSB_theta = Sigma_theta \ b_theta'; % D2 X 1 vector
end

if flagMH %Train's MH
    detSigma_theta = det(Sigma_theta);
    qJ = [];
    qNJ = [];
    qSum = [];
    qNorm = [];
    X1AvgQNorm = [];
    X2AvgQNorm = [];
    X1TimesX1 = [];
    X2TimesX2  = [];
else %Polya-Gamma MH
    detSigma_theta = [];
    if flagInit
        % number of visits to option J, for a certain (t, i)
        qJ = squeeze(Q(:, J, :)); % T X N matrix
        % number of visits to all the other options, for a certain (t, i)
        qNJ = squeeze(sum(Q(:, 1 : (J - 1), :), 2));
        % add together
        qSum = qJ + qNJ;
        % probability of each class given that it is not J
        qNorm = zeros(T, J - 1, N);
        for i = 1 : N
            for t = 1 : T
                qNorm(t, :, i) = Q(t, 1 : (J - 1), i) / (qNJ(t, i) + realmin);
            end
        end
        % X averaged on qNorm
        X1AvgQNorm = zeros(T, D1, N);
        X2AvgQNorm = zeros(T, D2, N);
        for i = 1 : N
            for t = 1 : T
                idr = ((i - 1) * TJ + (t - 1) *  J + 1) : ((i - 1) * TJ + t * J - 1); % 1, ..., J - 1
                idc = (i - 1) * D1 + 1 : i * D1;
                X1AvgQNorm(t, :, i) = qNorm(t, :, i) * X1S(idr, idc);
                X2AvgQNorm(t, :, i) = qNorm(t, :, i) * X2(1 : J - 1, :); % 1 X D2 vector
            end
        end
        % product of X
        X1TimesX1 = zeros(D1, D1, T, N);
        X2TimesX2 = zeros(D2, D2, T, N);
        for i = 1 : N
            for t = 1 : T
                X1TimesX1(:, :, t, i) = X1AvgQNorm(t, :, i)' * X1AvgQNorm(t, :, i);
                X2TimesX2(:, :, t, i) = X2AvgQNorm(t, :, i)' * X2AvgQNorm(t, :, i);
            end
        end
    else
        qJ = pcptPara.qJ;
        qNJ = pcptPara.qNJ;
        qSum = pcptPara.qSum;
        qNorm = pcptPara.qNorm;
        X1AvgQNorm = pcptPara.X1AvgQNorm;
        X2AvgQNorm = pcptPara.X2AvgQNorm;
        X1TimesX1 = pcptPara.X1TimesX1;
        X2TimesX2  = pcptPara.X2TimesX2;
    end
end

if flagInit == 1
    prod1 = [];
    prod2 = [];
    MNLLL = [];
else
    [prod1, prod2] = PrecProd(X1S, X2, modelPara.beta, modelPara.theta, dimPara.N, dimPara.T, dimPara.J);
    if flagMH == 1
        MNLLL = MNLLikelihood(prod1, prod2, Q, [], [], 1, 3);
    else
        MNLLL = [];
    end
end

pcptPara = PackPcptPara(...
    detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
    detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL);

return;