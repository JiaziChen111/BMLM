clear; clc; close all;
% rng('default');

%% Set parameters
flagDebug = 1;

% Data files
dataRoot = fullfile(pwd, 'Data');
dataFiles = dir(dataRoot);
numData = 0;
for i = 1 : length(dataFiles)
    if isempty(strfind(dataFiles(i).name, '.'))
        numData = numData + 1;
        dataFolder{numData} = fullfile(dataRoot, dataFiles(i).name);
        outPath{1, numData} = ['Output_', dataFiles(i).name];
    end
end
clear dataFiles;

flagMH = repmat([0, 0, 1, 0], 1, numData); %1 = Train, 0 = Polya-Gamma
flagSP = repmat([0, 0, -1, -1], 1, numData); % -1 = normal prior on theta, 0 = sparse prior on theta, 1 = sparse prior on b_theta
flagCS = repmat([1, 0, 0, 0], 1, numData); % 1 = common sparse prior, 0 = customer-specific sparse prior. this switch is ignored if flagSP = -1.
numRun = length(flagMH); % totoal number of runs
numExp = numRun / numData; % number of exprements for each data file
outPath = reshape(repmat(outPath, numExp, 1), [], 1)';

rngSeed = [];
n_burnin = 0;
n_MH = 1;

if flagDebug
    runs = [1];
    n_collect = 2e1 * ones(1, numRun); % samples to be saved
    iterDisplay = 1;
else
    runs = [1: numRun];
    n_collect = 2e3* ones(1, numRun); % samples to be saved%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    iterDisplay = 100;
end



%% Generate data and initialize beta
for i = 1 : numData
    for j = 1 : numExp
        iRun = (i - 1) * numExp + j;
        if ~isempty(strfind(dataFolder{i}, 'Toy'))
            filePrefix = 'PizzaData_N100_J1000';
        else
            filePrefix = 'Data_TX_08_10_AllStores_FrozenVeg_Final'; %'Sample';%['PizzaData_N', num2str(N(iRun)), '_J', num2str(J(iRun))];
        end
        dataFiles{iRun} = fullfile(dataFolder{i}, [filePrefix, '.mat']);
        initFiles{iRun} = [];%fullfile(dataFolder{i}, ['ParaInitValue_', filePrefix, '.csv']);;%['ParaInitValue_', filePrefix, '.csv'];
        if ~isempty(initFiles{iRun}) && ~exist(initFiles{iRun}, 'file')
            FitIndividualFixLogit(dataFiles{iRun});
        end
    end
end



%% Fit MNL-DPM
% From now on all parallel computing in the called functions should be disabled.
flagPar = 0;
if flagDebug == 0 
    delete(gcp);
    parpool(min(2, length(runs))); % max number of threads (depending on memory usage)
    ms.UseParallel = true;
end

if flagDebug
    for iRun = runs
%         fprintf(sprintf('Run %d: J = %d, N = %d, flagMH = %d, flagSP = %d, flagCS = %d\n', iRun, J(iRun), N(iRun), flagMH(iRun), flagSP(iRun), flagCS(iRun)));
        optPara = PackOptPara(dataFiles{iRun}, initFiles{iRun}, outPath{iRun}, flagMH(iRun), flagSP(iRun), flagCS(iRun), flagPar, rngSeed, n_burnin, n_collect(iRun), n_MH, iterDisplay);
        DPM_MNL(optPara);
    end
else
    fid = fopen('log.txt', 'a');
    %     s = ['\nStart running ', strrep(pwd, '\', '\\'), '\t', datestr(datetime), '\n'];
    s = ['\nStart running ', strrep(pwd, '\', '\\'), '\t', datestr(clock), '\n'];
    fprintf(fid, s);
    fclose(fid);
    parfor i = 1 : length(runs)
        iRun = runs(i);
        optPara = PackOptPara(dataFiles{iRun}, initFiles{iRun}, outPath{iRun}, flagMH(iRun), flagSP(iRun), flagCS(iRun), flagPar, rngSeed, n_burnin, n_collect(iRun), n_MH, iterDisplay);
        fullOutPath = fullfile(outPath{iRun}, DPM_MNL(optPara));
        fid = fopen('log.txt', 'a');
        %         s = [num2str(iRun), ': Finished ', strrep(dataFiles{iRun}, '\', '\\'), ', saved in ', strrep(fullOutPath, '\', '\\'), '\t', datestr(datetime), '\n'];
        s = [num2str(iRun), ': Finished ', strrep(dataFiles{iRun}, '\', '\\'), ', saved in ', strrep(fullOutPath, '\', '\\'), '\t', datestr(clock), '\n'];
        fprintf(fid, s);
        fclose(fid);
    end
end

if flagDebug == 0
    delete(gcp);
end

return;
