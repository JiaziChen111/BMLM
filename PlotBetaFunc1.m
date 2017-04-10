function [] = PlotBetaFunc1( flagMH, flagSP,flagCS, dataset  )
%Function to plot beta's for every dataset

J = 1e3;
N = 1e2;
% flagMH = 0;
% flagCS = 1;
n_burnin = 0;
n_collect = 2e3;
% n_collect = 200;

n_MH = 1;
flagCoda = 0;
% Thinning factor
thin = 1 %4;
ids = 1001 %3001;
% ids = 1 %3001;

rootPath = pwd;
% dataset = 'ToyCustomerSparsity09'; % 'ToyNoSparsity';%
dataFile = ['PizzaData_', sprintf('N%d_J%d', N, J), '.csv'];

%% Load true beta
A = load(fullfile('Data', dataset, [dataFile(1 : end-4), '_all.mat']));

%% Load estimate of the new model
% flagSP = 0;
% Result folder
outPath = fullfile(rootPath, ['Output_',dataset]); 
optPara = PackOptPara(dataFile, [], outPath, flagMH, flagSP, flagCS, [], [], n_burnin, n_collect, n_MH, []);
outPath = GenOutpath(optPara);
outPath1 = outPath;
% Load
load(fullfile(outPath, 'Samples_beta.mat'));
betaEst1 = mean(beta_collection(:, :, ids : thin : n_collect), 3);

%% Load estimate of the old model
% flagSP = -1;
% % Result folder
% outPath = fullfile(rootPath, ['Output_',dataset]); 
% optPara = PackOptPara(dataFile, [], outPath, flagMH, flagSP, flagCS, [], [], n_burnin, n_collect, n_MH, []);
% outPath = GenOutpath(optPara);
% % Load
% load(fullfile(outPath, 'Samples_beta.mat'));
% betaEst2 = mean(beta_collection(:, :, ids : thin : n_collect), 3);

%% Plot density
D = size(A.beta, 2);

figure(100),
for j = 1 : D
    subplot(2, ceil(D/2), j)
    [f, xi] = ksdensity(A.beta(: ,j), ...
        'function', 'pdf', 'width', 1);
    plot(xi, f, 'b-');
    hold on
    [f, xi] = ksdensity(betaEst1(: ,j), ...
        'function', 'pdf', 'width', 1);
    plot(xi, f, 'r--');
%     [f, xi] = ksdensity(betaEst2(: ,j), ...
%         'function', 'pdf', 'width', 1);
%     plot(xi, f, 'k:');
    xlabel(['\beta_', num2str(j)]);
    ylabel('density');
    hold off
   
end
saveas(gca,[outPath1,filesep,'Beta'],'epsc');

end

