function [detSigma_beta, invSigma_beta, LSigma_beta, invSB_beta, NLLBeta, invSigma_theta, invSB_theta, invTau_theta, ...
    detSigma_theta, qJ, qNJ, qSum, qNorm, X1AvgQNorm, X2AvgQNorm, X1TimesX1, X2TimesX2, prod1, prod2, MNLLL] = UnpackPcptPara(pcptPara)

detSigma_beta = pcptPara.detSigma_beta; 
invSigma_beta = pcptPara.invSigma_beta; 
LSigma_beta = pcptPara.LSigma_beta; 
invSB_beta = pcptPara.invSB_beta;
NLLBeta = pcptPara.NLLBeta;
invSigma_theta = pcptPara.invSigma_theta;
invSB_theta = pcptPara.invSB_theta;
invTau_theta = pcptPara.invTau_theta;
detSigma_theta = pcptPara.detSigma_theta;
qJ = pcptPara.qJ;
qNJ = pcptPara.qNJ; 
qSum = pcptPara.qSum; 
qNorm = pcptPara.qNorm; 
X1AvgQNorm = pcptPara.X1AvgQNorm; 
X2AvgQNorm = pcptPara.X2AvgQNorm; 
X1TimesX1 = pcptPara.X1TimesX1; 
X2TimesX2 = pcptPara.X2TimesX2;
prod1 = pcptPara.prod1;
prod2 = pcptPara.prod2;
MNLLL = pcptPara.MNLLL;

return;