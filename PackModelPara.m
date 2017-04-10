function modelPara = PackModelPara(c, C, beta, b_beta, Sigma_beta, theta, b_theta, Sigma_theta, tau_theta, lambda_theta)

modelPara.c = c;
modelPara.C = C ;

modelPara.beta = beta;
modelPara.b_beta = b_beta;
modelPara.Sigma_beta = Sigma_beta;

modelPara.theta = theta;
modelPara.b_theta = b_theta;
modelPara.Sigma_theta = Sigma_theta;
modelPara.tau_theta = tau_theta;
modelPara.lambda_theta = lambda_theta;

return;