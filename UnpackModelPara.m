function [c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta] = UnpackModelPara(modelPara)

c = modelPara.c;
C = modelPara.C ;

beta = modelPara.beta;
b_beta = modelPara.b_beta;
Sigma_beta = modelPara.Sigma_beta;

theta = modelPara.theta;
b_theta = modelPara.b_theta;
Sigma_theta = modelPara.Sigma_theta;
tau_theta = modelPara.tau_theta;
lambda_theta = modelPara.lambda_theta;

return;