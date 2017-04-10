function hyperPara = SetHyperPara(dimPara, flagSP)

if nargin < 2
    flagNIW = 1;
end

[N, T, J, D1, D2, NT, TJ, NTJ] = UnpackDimPara(dimPara);

% beta: DPM hyper-parameters
alpha = 10;

if (1)
    b_0beta = zeros(1, D1);
    upsilon_0beta = D1 + 2;
    S_0beta = eye(D1);% * (upsilon_0beta - D1 - 1); % Note: Matlab's IW is different from that in Train (2003). Need to adjust parameters accordingly.
    
    LS_0beta = inv(chol(S_0beta))';
    kappa_0beta = 1;%0.1;
    Sigma_0beta = 10 * eye(D1); %S_0beta / kappa_0beta;
    invSigma_0beta = inv(Sigma_0beta);
    LSigma_0beta = chol(Sigma_0beta, 'lower');
else
    % fix covariance Sigma_theta in each cluster
    % b_beta prior N(b_0beta, S_0beta)
    b_0beta = zeros(1, D1);
    upsilon_0beta = [];
    S_0beta = 10 * eye(D1);% * (upsilon_0beta - D1 - 1); % Note: Matlab's IW is different from that in Train (2003). Need to adjust parameters accordingly.   
    LS_0beta = inv(chol(S_0beta))';
    kappa_0beta = [];
    
    Sigma_0beta = eye(D1);
    invSigma_0beta = inv(Sigma_0beta);
    LSigma_0beta = chol(Sigma_0beta, 'lower');
end

% theta: normal hyper-parameters
if flagSP == 1 %sparse prior on b_theta
    upsilon_0theta = D2 + 2;
    kappa_0theta = 0.01;
    b_0theta = [];
    S_0theta = [];
    % for lapalacian prior
    r_theta0 = 1;
    delta_theta0 = 0.1;
end

if flagSP == 0 %sparse prior on theta so we don't need these two hyper-parameters
    upsilon_0theta = [];
    kappa_0theta = [];
    b_0theta = [];
    S_0theta = [];
    % for lapalacian prior
    r_theta0 = 1;
    delta_theta0 = 0.1;
end


if flagSP == -1
    % theta: normal hyper-parameters
    b_0theta = zeros(1.0, D2);
    upsilon_0theta = D2 + 2;
    S_0theta = eye(D2);% * (upsilon_0theta - D2 - 1); % Note: Matlab's IW is different from that in Train (2003). Need to adjust parameters accordingly.
    %     LS_0theta = inv(chol(S_0theta))';
    kappa_0theta = 0.1;%0.1;
    %     Sigma_0theta = S_0theta / (upsilon_0theta - D2 - 1) / kappa_0theta;
    %     invSigma_0theta = []; % avoid computing the inverse of a big matrix
    %     LSigma_0theta = chol(Sigma_0theta, 'lower');
    r_theta0 = [];
    delta_theta0 = [];
end

% Pack hyper-parameters
hyperPara = PackHyperPara(alpha, ...
    b_0beta, kappa_0beta, upsilon_0beta, S_0beta, LS_0beta, Sigma_0beta, invSigma_0beta, LSigma_0beta, ...
    b_0theta, kappa_0theta, upsilon_0theta, S_0theta, r_theta0, delta_theta0);

return;