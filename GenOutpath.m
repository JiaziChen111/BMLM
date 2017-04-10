function outPath = GenOutpath(optPara)

[dataFile, initFile, outPath, flagMH, flagSP, flagCS, flagPar, rngSeed, n_burnin, n_collect, n_MH, prtIntv] = UnpackOptPara(optPara);

if isempty(outPath)
    rootPath = 'Output';
else
    rootPath = outPath;
end
% dataFile could be file name only or full path
[pathstr, name, ext] = fileparts(dataFile);
outPath = name(1 : min(23, length(name)));
if flagMH == 1
    outPath = [outPath, '_Train'];
else
    outPath = [outPath, '_PG'];
end
if flagSP > -1
    if flagCS == 1
        outPath = [outPath, '_Common'];
    else
        outPath = [outPath, '_Indv'];
    end
end
switch flagSP 
    case 1
        outPath = [outPath, '_SparseOn_b_theta'];
    case 0
        outPath = [outPath, '_SparseOn_theta'];
    otherwise
        outPath = [outPath, '_NoSparse'];
end
outPath = sprintf([outPath, '_nBurnin_%d_nCollect_%d_nMH_%d'], n_burnin, n_collect, n_MH);
outPath = fullfile(rootPath, outPath);
return;