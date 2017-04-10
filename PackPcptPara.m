function pcptPara = PackPcptPara(detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
    detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL)

pcptPara.detSigma_beta = detSigma_beta; 
pcptPara.invSigma_beta = invSigma_beta; 
pcptPara.LSigma_beta = LSigma_beta; 
pcptPara.invSB_beta = invSB_beta;
pcptPara.NLLBeta = NLLBeta;
pcptPara.invSigma_theta = invSigma_theta;
pcptPara.invSB_theta = invSB_theta;
pcptPara.invTau_theta = invTau_theta;
pcptPara.detSigma_theta = detSigma_theta;
pcptPara.qJ = qJ;
pcptPara.qNJ = qNJ; 
pcptPara.qSum = qSum; 
pcptPara.qNorm = qNorm; 
pcptPara.X1AvgQNorm = X1AvgQNorm; 
pcptPara.X2AvgQNorm = X2AvgQNorm; 
pcptPara.X1TimesX1 = X1TimesX1; 
pcptPara.X2TimesX2 = X2TimesX2;
pcptPara.prod1 = prod1;
pcptPara.prod2 = prod2;
pcptPara.MNLLL = MNLLL;

return;