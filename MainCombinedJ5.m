%% ToySparsity00
% 
%% PG Indiv Sparse on Theta

MainPostFunc(5,100,0,0,0,'ToySparsity00');
PlotBetaFunc(5,100,00,0,'ToySparsity00');
PlotThetaFunc(5,100,00,0,'ToySparsity00');

%% PG Common Sparse on Theta

MainPostFunc(5,100,00,1,'ToySparsity00');
PlotBetaFunc(5,100,00,1,'ToySparsity00');
PlotThetaFunc(5,100,00,1,'ToySparsity00');

%% Train No Sparse

MainPostFunc(5,100,1,-1,0,'ToySparsity00');
PlotBetaFunc(5,100,1,-1,0,'ToySparsity00');
PlotThetaFunc(5,100,1,-1,0,'ToySparsity00');

%% PG No Sparse

MainPostFunc(5,100,0,-1,0,'ToySparsity00');
PlotBetaFunc(5,100,0,-1,0,'ToySparsity00');
PlotThetaFunc(5,100,0,-1,0,'ToySparsity00');


