function errCode = CheckMembership(C, c, clusSize)
% check if the component size is consistent with the membership variables
% errCode = 0: no error;

errCode = 0;

if C ~= length(unique(c))
    errCode = 100;
    return;
end

if C ~= length(clusSize)
    errCode = 200;
    return;
end

if length(clusSize) ~= length(unique(c))
    errCode = 300;
    return;
end

for k = 1 : C
    if sum(c==k) ~= clusSize(k)
        errCode = k;
        return;
    end
end



