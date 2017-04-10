function optPara = PackOptPara(dataFile, initFile, outPath, flagMH, flagSP, flagCS, flagPar, rngSeed, n_burnin, n_collect, n_MH, prtIntv)

% Data file name
optPara.dataFile = dataFile;
% Parameter initilization
optPara.initFile = initFile;
% Where to save the results
optPara.outPath = outPath;
% Which MH method to use: 1 = Train, 0 = Polya-Gamma
optPara.flagMH = flagMH;
% Where to put sparse prior: 1 = b_theta, 0 = theta
optPara.flagSP = flagSP;
% Common or individulized sparse prior: 1 = common, 0 = individualized
optPara.flagCS = flagCS;
% Parallel computing
optPara.flagPar = flagPar;
% Seed for random number generator
optPara.rngSeed = rngSeed;
% Gibbs sampling parameters
optPara.n_burnin = n_burnin;
optPara.n_collect = n_collect; 
% Number of MH steps in each Gibbs iteration - let's keep it as 1 for now
optPara.n_MH = n_MH;
% Print interval
optPara.prtIntv = prtIntv;

return;