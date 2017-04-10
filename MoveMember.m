function [compPara, c] = MoveMember(compPara, i, l, c)

% Move member i from the current component to component l
k = c(i);
compPara.compSize(k) = compPara.compSize(k) - 1;

c(i) = l;
compPara.compSize(l) = compPara.compSize(l) + 1;

% Drop component k if now it's empty
if compPara.compSize(k) == 0
    [compPara, c] = DropComponent(compPara, k, c);
end

return
