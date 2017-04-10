function [alpha, b_0beta, kappa_0beta, upsilon_0beta, S_0beta, LS_0beta, Sigma_0beta, invSigma_0beta, LSigma_0beta, ...
    b_0theta, kappa_0theta, upsilon_0theta, S_0theta, r_0theta, delta_0theta] ...
    = UnpackHyperPara(hyperPara)

% beta: DPM hyper-parameters
alpha = hyperPara.alpha; 
b_0beta = hyperPara.b_0beta; 
kappa_0beta = hyperPara.kappa_0beta; 
upsilon_0beta = hyperPara.upsilon_0beta; 
S_0beta = hyperPara.S_0beta; 
LS_0beta = hyperPara.LS_0beta; 
Sigma_0beta = hyperPara.Sigma_0beta; 
invSigma_0beta = hyperPara.invSigma_0beta; 
LSigma_0beta = hyperPara.LSigma_0beta; 

% theta: normal hyper-parameters
b_0theta = hyperPara.b_0theta; 
kappa_0theta = hyperPara.kappa_0theta; 
upsilon_0theta = hyperPara.upsilon_0theta; 
S_0theta = hyperPara.S_0theta; 
% LS_0theta = hyperPara.LS_0theta; 
% Sigma_0theta = hyperPara.Sigma_0theta; 
% invSigma_0theta = hyperPara.invSigma_0theta; 
% LSigma_0theta = hyperPara.LSigma_0theta; 

r_0theta = hyperPara.r_0theta; 
delta_0theta = hyperPara.delta_0theta;

return;