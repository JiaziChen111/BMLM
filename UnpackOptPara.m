function [dataFile, initFile, outPath, flagMH, flagSP, flagCS, flagPar, rngSeed, n_burnin, n_collect, n_MH, prtIntv] = UnpackOptPara(optPara)

% Data file name
dataFile = optPara.dataFile;
% Parameter initilization
initFile = optPara.initFile;
% Where to save the results
outPath = optPara.outPath;
% Which MH method to use: 1 = optPara.Train, 0 = optPara.Polya-Gamma
flagMH = optPara.flagMH;
% Where to put sparse prior: 1 = optPara.b_theta, 0 = optPara.theta
flagSP = optPara.flagSP;
% Common or individulized sparse prior: 1 = common, 0 = individualized
flagCS = optPara.flagCS;
% Parallel computing
flagPar = optPara.flagPar;
% Seed for random number generator
rngSeed = optPara.rngSeed;
% Gibbs sampling parameters
n_burnin = optPara.n_burnin;
n_collect = optPara.n_collect; 
% Number of MH steps in each Gibbs iteration - let's keep it as 1 for now
n_MH = optPara.n_MH;
% Print interval
prtIntv = optPara.prtIntv;

return;