function [compPara, cNew] = DropComponent(compPara, k, c)

% Drop the kth component (here it's assumed to be an empty component)
% Adjust index for the members in component [k+1, C]

compPara.compSize(k) = [];
compPara.b_beta(k, :) = [];

if ~isempty(compPara.Sigma_beta)
    compPara.Sigma_beta(:, :, k) = [];
end

if ~isempty(compPara.invSigma_beta)
    compPara.invSigma_beta(:, :, k) = [];
end

if ~isempty(compPara.detSigma_beta)
    compPara.detSigma_beta(k) = [];
end

cNew = c;
for l = k + 1 : compPara.C
    cNew(find(c == l)) = l - 1; 
end

compPara.C = compPara.C - 1;

return
