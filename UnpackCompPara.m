function [C, compSize, b_beta, Sigma_beta, invSigma_beta, detSigma_beta] = UnpackCompPara(compPara)

C = compPara.C;
compSize = compPara.compSize;
b_beta = compPara.b_beta;
Sigma_beta = compPara.Sigma_beta;
invSigma_beta = compPara.invSigma_beta;
detSigma_beta = compPara.detSigma_beta; 

return;