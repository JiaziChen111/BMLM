function x = MultVMV(v1, M, v2, flagInv)
% Efficient way of computing the product of a vector, a matrix and a second vector 
% Setting 1:
% v1: N X D vector
% M:  D X D X N matrix
% v2: N X D vector
% x:  N X 1 vector. x(i) = v1(i, :) * M(:, :, i) * v2(i, :)'
% Setting 2:
% v1: N X D vector
% M:  D X D X N matrix
% x:  N X D vector. x(i, :) = v1(i, :) * M(:, :, i)'
% flagInv = 1: v1 * inv(M) * v2

if nargin < 4
    flagInv = 0;
end

if nargin > 2
    flagV2 = 1;
else
    flagV2 = 0;
end

[N, D] = size(v1);
dimM = size(M);

if dimM(1) ~= D || dimM(2) ~= D
    fprintf('Error: Check the size of input arguments!');
end

if length(dimM) == 3
    if dimM(3) ~= N
        fprintf('Error: Check the size of input arguments!');
    end
end

if flagV2
    dimV2 = size(v2);
    if dimV2(1) ~= N || dimV2(2) ~= D
        fprintf('Error: Check the size of input arguments!');
    end
end

% v1 is a scalar
if N == 1 && D == 1
    if flagV2
        if flagInv
            x = v1 / M * v2;
        else
            x = v1 * M * v2;
        end
    else
        if flagInv
            x = v1 / M;
        else
            x = v1 * M;
        end
    end
end

% v1 is a row vector
if N == 1 && D > 1
    if flagV2
        if flagInv
            x = v1 * (M \ v2');
        else
            x = v1 * M * v2';
        end
    else
        if flagInv
            x = v1 * inv(M);
        else
            x = v1 * M;
        end
    end
end

% v1 is a column vector
if N > 1 && D == 1
    if flagV2
        if flagInv
            x = v1 ./ squeeze(M) .* v2;
        else
            x = v1 .* squeeze(M) .* v2;
        end
    else
        if flagInv
            x = v1 ./ squeeze(M);
        else
            x = v1 .* squeeze(M);
        end
    end
end

% v1 is a matrix
if N > 1 && D > 1
    if flagInv && flagV2
        for i = 1 : N
             x(i, :) = v1(i, :) * (M(:, :, i) \ v2(i, :)');
        end
    else
        v1Rep = repmat(reshape(v1', [D, 1, N]), [1, D, 1]); % D X D X N
        x = sum(v1Rep .* M, 1); % 1 X D X N
        if flagV2
            x = sum(squeeze(x)' .* v2, 2); % N X 1
        else
            x = squeeze(x)';
        end
    end
end

return;
