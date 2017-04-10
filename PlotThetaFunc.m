function [ ] = PlotThetaFunc( J,N,flagMH, flagSP,flagCS, dataset )
%Function to plot theta's for every dataset

% J = 1e3;
% N = 1e2;
% flagMH = 1;
% flagSP = -1;
% flagCS = 0;
n_burnin = 0;
n_collect = 2e3;
n_MH = 1;
flagCoda = 0;
% Thinning factor
thin = 1 %4;
ids = 1001; %3001;
rootPath = pwd;
dataFile = ['PizzaData_', sprintf('N%d_J%d', N, J), '.csv'];
% dataset = 'ToyCustomerSparsity09'; %'ToyNoSparsity'; %

%% Load true theta
load(fullfile('Data', dataset, [dataFile(1 : end-4), '_all.mat']));
meanTheta = mean(theta, 1);

%% Load estimate of the new model
% flagSP = 0;
% Result folder

outPath = fullfile(rootPath, ['Output_',dataset]);
optPara = PackOptPara(dataFile, [], outPath, flagMH, flagSP, flagCS, [], [], n_burnin, n_collect, n_MH, []);
outPath = GenOutpath(optPara);
outPath1 = outPath;
% Load
load(fullfile(outPath, 'Samples_theta.mat'));
meanThetaEst1 = squeeze(mean(mean(theta_collection(:, :, 11:20), 1), 3));%ids : thin : n_collect), 1), 3));
err1 = sqrt(mean(mean((theta_collection(:, :, 11:20) - repmat(theta, [1 1 10])).^2, 3), 1));

%% Load estimate of the old model
flagSP = -1;
% Result folder
outPath = fullfile(rootPath, ['Output_',dataset]);
optPara = PackOptPara(dataFile, [], outPath, flagMH, flagSP, flagCS, [], [], n_burnin, n_collect, n_MH, []);
outPath = GenOutpath(optPara);
% Load
load(fullfile(outPath, 'Samples_theta.mat'));
meanThetaEst2 = squeeze(mean(mean(theta_collection(:, :, 11:20), 1), 3));%ids : thin : n_collect), 1), 3));
err2 = sqrt(mean(mean((theta_collection(:, :, 11:20) - repmat(theta, [1 1 10])).^2, 3), 1));

%% Plot density
figure,
subplot(3, 1, 1)
plot(meanTheta, 'b-');
subplot(3, 1, 2)
plot(meanThetaEst1, 'r--');
subplot(3, 1, 3)
plot(meanThetaEst2, 'k:');
saveas(gca,[outPath1,filesep,'Theta'],'epsc');

%% Load estimate of the new model
% flagSP = 0;
% % Result folder
% outPath = fullfile(rootPath, ['Output_',dataset]);
% optPara = PackOptPara(dataFile, [], outPath, flagMH, flagSP, flagCS, [], [], n_burnin, n_collect, n_MH, []);
% outPath = GenOutpath(optPara);
% % Load
% load(fullfile(outPath, 'Samples_theta.mat'));
% err1 = sqrt(mean(mean((theta_collection(:, :, 11:20) - repmat(theta, [1 1 10])).^2, 3), 1));

%% Load estimate of the old model
% flagSP = -1;
% % Result folder
% outPath = fullfile(rootPath, ['Output_',dataset]);
% optPara = PackOptPara(dataFile, [], outPath, flagMH, flagSP, flagCS, [], [], n_burnin, n_collect, n_MH, []);
% outPath = GenOutpath(optPara);
% % Load
% load(fullfile(outPath, 'Samples_theta.mat'));
% err2 = sqrt(mean(mean((theta_collection(:, :, 11:20) - repmat(theta, [1 1 10])).^2, 3), 1));
% 
%% Plot density
figure,
subplot(2, 1, 1)
plot(err1, 'r--');
subplot(2, 1, 2)
plot(err2, 'k:');
saveas(gca,[outPath1,filesep,'thErr'],'epsc');


mean(err1)
mean(err2)



end

