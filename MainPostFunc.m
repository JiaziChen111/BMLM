function[]  = MainPostFunc( J,N,flagMH, flagSP,flagCS, dataset )
% A function to perform the same operations as MainPost.m

% J = 1e3
% N = 1e2;
% flagMH = 1;
% flagSP = -1;
% flagCS = 0;
n_burnin = 0;
n_collect = 2e3;
% n_collect = 200;

n_MH = 1;
flagCoda = 0;
% Thinning factor
thin = 1; %4;
ids = n_collect / 2 + 1; %3001;

% dataset = 'ToyCustomerSparsity09'; %'ToyNoSparsity'; %
dataFile = ['PizzaData_', sprintf('N%d_J%d', N, J)]; %'.csv'%


% Result folder
%% Specify the path to the result folder
rootPath = pwd;
outPath = fullfile(rootPath, ['Output_',dataset]); 
optPara = PackOptPara(dataFile, [], outPath, flagMH, flagSP, flagCS, [], [], n_burnin, n_collect, n_MH, []);
outPath = GenOutpath(optPara);

%% Class allocation analysis
% First load original data
A = load(fullfile(rootPath,'Data', dataset, [dataFile, '_all.mat'])); %dataFile(1 : end-4)
if ~exist('c', 'var') % didn't save c in earlier versions of GeneratePizzaData()
    C = 18;
    D1 = size(A.beta, 2);
    % mean coefficients
    b_A.beta = zeros(C, D1); % not senstive to brand
    % mean coef for price
    b_A.beta(:, 1) = [ones(1, 6) * (-4), zeros(1, 6), ones(1, 6) * 4]';
    % for type
    b_A.beta(:, 2) = repmat([-5, -5, 0, 0, 5, 5],[1, 3]);
    % for walking distance to the pizza
    b_A.beta(:, 3) = repmat([-6, 0], [1, 9])';
    
    % covariance
    sigma = ones(C, D1) * 0.2;
    sigma(:, 1) = [ones(1, 6) * 0.5, ones(1, 6) * 1.5, ones(1, 6) * 0.5];
    sigma(:, 2) = repmat([1, 1, 2, 2, 1, 1],[1, 3]);
    sigma(:, 3) = repmat([0.5, 1], [1, 9])';
    Sigma_beta = zeros(D1, D1, C);
    for k = 1 : C
        Sigma_beta(:, :, k) = diag(sigma(k, :));
    end
    % weights of each latent class
    % marginal weights
    mw = {};
    mw{1} = [0.3, 0.6, 0.1];
    mw{2} = [0.5, 0.0, 0.5]; %[0.4, 0.2, 0.4];
    mw{3} = [0.2, 0.8];
    % weights
    w = zeros(C, 1);
    iter = 1;
    for i = 1 : 3
        for j = 1: 3
            for k = 1: 2
                w(iter) = mw{1}(i) * mw{2}(j) * mw{3}(k);
                iter = iter + 1;
            end
        end
    end
    %sum(w)
    % Find the nearest latent class to the true A.beta
    c =  AssignClass(A.beta, A.b_beta);
else
    C = 18;
end

% Change name
CTrue = C;
cTrue = A.c;
b_betaTrue = A.b_beta;
betaTrue = A.beta;
thetaTrue = A.theta;
% 
figure,
hist(cTrue, 1 : CTrue);
title('Latent class allocation - true');
saveas(gca,[outPath,filesep,'LCTrue'],'epsc');

% Show class center
[(1:18)', b_betaTrue(:, 1 : 3)]

% Load results
load(fullfile(outPath, 'Samples_c.mat'));
load(fullfile(outPath, 'Samples_beta.mat'));
load(fullfile(outPath, 'Samples_b_beta.mat'));
load(fullfile(outPath, 'Samples_Sigma_beta.mat'));
load(fullfile(outPath, 'Samples_theta.mat'));
load(fullfile(outPath, 'Samples_b_theta.mat'));
load(fullfile(outPath, 'Samples_Sigma_theta.mat'));
load(fullfile(outPath, 'Samples_lambda_theta.mat'));
load(fullfile(outPath, 'Samples_tau_theta.mat'));

% Take the last sample
c = c_collection(:, end);
b_beta = b_beta_collection{end};
Sigma_beta = Sigma_beta_collection{end};
C = size(b_beta, 1);

% figure,
% hist(c, 1 : C);
% title('Latent class allocation - est');
% saveas(gcf,[outPath,filesep, 'LCEst'],'epsc');


%% Convergence analysis

% number of clusters
n_collect = size(c_collection, 2);
for iter = 1:n_collect
    C(1, iter) = length(unique(c_collection(:,iter)));
end
flagKDE = 0;
[mC, prRL_C, prG_C] = EvalSamp(C, thin, ids, [], flagCoda, flagKDE, 'number of clusters');

% A.beta
for i = 1 : 10%N
    figure,
    plot(squeeze(beta_collection(i,1,:)));
    xlabel(sprintf('\beta_{%d,1}', i));
    saveas(gcf,[outPath,filesep, sprintf('b_%d', i)],'epsc');

end
flagKDE = 1;
[mbeta, prRL_beta, prG_beta] = EvalSamp(beta_collection, thin, ids, [], flagCoda, flagKDE, 'beta'); % do KDE for beta

% tau_theta
figure,
plot(mean(tau_theta_collection(:,:,end), 1));
saveas(gcf,[outPath,filesep, 'tau'],'epsc');
flagKDE = 0;
[mTheta, prRL_b_theta, prG_b_theta] = EvalSamp(tau_theta_collection, thin, ids, [], flagCoda, flagKDE, 'tau_theta');

% Plot theta related
if flagSP == 1
    figure(201)
    subplot(3, 1, 1)
    plot(theta_collection(1,:,end));
    ylabel('\theta_{i=1}');
    subplot(3, 1, 2)
    plot(b_theta_collection(1, :, end));
    ylabel('b_{\theta}');
    subplot(3, 1, 3)
    Sigma_theta = Sigma_theta_collection(:, :, end);
    plot(diag(Sigma_theta));
    ylabel('diag(\Sigma_{\theta})');
else
    % convert to binary matrix
    randomset = tau_theta_collection(:, :, end);
    figure(201)
    imagesc(randomset);
end
saveas(gcf,[outPath,filesep,'tplot'],'epsc');

% computer error
betaEst = mean(beta_collection(:, :, ids:thin:end), 3);
err = mean(sqrt(sum((betaEst - betaTrue).^2, 2)));
fprintf(['error = ', num2str(err), '\n']);


% check on theta
thetaEst = theta_collection(:, :, end);
load(fullfile('Data', dataset, [dataFile, '.mat'])); %dataFile(1 : end-4)
figure, 
subplot(2, 1, 1)
plot(thetaTrue(2, :) .* Z(2, 1 : end-1));
subplot(2, 1, 2)
plot(thetaEst(2, :) .* Z(2, 1 : end-1));

saveas(gcf,[outPath,filesep,'tCplot'],'epsc');

meanTheta = zeros(1, D2);
for i = 1 : N
    meanTheta = meanTheta + thetaEst(i, :) .* Z(i, 1 : end-1);
end
meanTheta = meanTheta / N;
figure,
plot(meanTheta);
saveas(gcf,[outPath,filesep, 'tMplot'],'epsc');





end

