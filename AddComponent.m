function compPara = AddComponent(compPara, b_betaNew, Sigma_betaNew)
% Add an empty new component


compPara.C = compPara.C + 1;
compPara.compSize = [compPara.compSize; 0];
compPara.b_beta = [compPara.b_beta; b_betaNew];
compPara.Sigma_beta = cat(3, compPara.Sigma_beta, Sigma_betaNew);
compPara.invSigma_beta = cat(3, compPara.invSigma_beta, inv(Sigma_betaNew));
compPara.detSigma_beta = [compPara.detSigma_beta; det(Sigma_betaNew)];
