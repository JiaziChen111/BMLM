function L = MNLLikelihood(X_1, X_2, Q, beta, theta, flagLog, flagFast)
% Compute the likelihood of a MNL model

% flagLog:  Scalar. 1 = in log scale, 0 = in normal scale.
% flagFast: Use recomputed X_1 * beta and/or X_2 * theta.

% If flagFast == 0: no pre-computing
% X_1:              (N*TJ) X (N*D1) sparse matrix 
% X_2:              J X D2 matrix
% Q:                T X J X N matrix
% beta:             N X D1 matrix
% theta:            N X D2 matrix

% If flagFast == 1: use pre-computed prod1
% X_1:              prod_1, J X T X N matrix
% X_2:              J X D2 matrix
% Q:                T X J X N matrix
% beta:             []
% theta:            N X D2 matrix

% If flagFast == 2: use pre-computed prod2
% X_1:              (N*TJ) X (N*D1) sparse matrix 
% X_2:              prod_2, J X T X N matrix
% Q:                T X J X N matrix
% beta:             N X D1 matrix
% theta:            []

% If flagFast == 3: use pre-computed prod1 and prod2
% X_1:              prod_1, J X T X N matrix
% X_2:              prod_2, J X T X N matrix
% Q:                T X J X N matrix
% beta:             []
% theta:            []


if nargin < 7
        fprintf('Error: All the input arguments have to be specified. Some of them could be empty though.')
end
    
if flagLog ~= 1 && flagLog ~= 0
    fprintf('Error: flagLog can only be 0 or 1!');
    return;
end


% dimX1 = dim(X1);
% if dimX1(1) ~= N || dimX1(2) ~= T || dimX1(3) ~= J 
%     fprint('Check X1!');
%     return;
% end

dims = size(Q);
T = dims(1);
J = dims(2);
if length(dims) == 3
    N  = dims(3);
else
    N = 1;
end

% slow loop
% logL = 0;
% for i = 1 : N
%     for t = 1 : T
%         a = exp( X_1( :, :, t, i) * beta( i, : )' + X_2 * theta( i, : )');
%         logL = logL + sum( Q( t, :, i ) * ( log(a) - log(sum(a)) ) );
%     end
% end


if N > 1
    switch flagFast
        case 0
            %prod1 = squeeze(sum(X_1 .* repmat(reshape(beta', [1, D1, 1, N]), [J, 1, T, 1]), 2)); % J X T X N
            prod1 = reshape(X_1 * reshape(beta', [], 1), [J, T, N]);
            prod2 = repmat(reshape([theta'; zeros(1, N)], [J, 1, N]), [1, T, 1]); % J X T X N
        case 1 % use pre-computed prod1
            prod1 = X_1;
            prod2 = repmat(reshape([theta'; zeros(1, N)], [J, 1, N]), [1, T, 1]); % J X T X N
        case 2 % use pre-computed prod2
            % prod1 = squeeze(sum(X_1 .* repmat(reshape(beta', [1, D1, 1, N]), [J, 1, T, 1]), 2)); % J X T X N
            prod1 = reshape(X_1 * reshape(beta', [], 1), [J, T, N]);
            prod2 = X_2;
        otherwise % use pre-computed prod1 and prod2
            prod1 = X_1;
            prod2 = X_2;
    end
    g = prod1 + prod2;
    logSumExpG = log(sum(exp(g), 1)); % 1 X T X N
    logL = permute(Q, [2, 1, 3]) .* (g - repmat(logSumExpG, [J, 1, 1])); % J X T X N
else
    switch flagFast
        case 0
            % prod1 = squeeze(sum(X_1 .* repmat(beta, [J, 1, T]), 2));
            prod1 = reshape(X_1 * beta, [J, T]);
            prod2 = repmat([theta'; 0], [1, T]);
        case 1 % use pre-computed prod1
            prod1 = X_1;
            prod2 = repmat([theta'; 0], [1, T]);
        case 2 % use pre-computed prod2
            % prod1 = squeeze(sum(X_1 .* repmat(beta, [J, 1, T]), 2));
            prod1 = reshape(X_1 * beta, [J, T]);
            prod2 = X_2;
        otherwise % use pre-computed prod1 and prod2
            prod1 = X_1;
            prod2 = X_2;
    end
    g = prod1 + prod2;
    logSumExpG = log(sum(exp(g), 1)); % 1 X T
    logL = Q' .* (g - repmat(logSumExpG, [J, 1])); % J X T
end

logL = squeeze(sum(sum(logL, 1), 2)); % N X 1

if flagLog == 0
    L = exp( logL );
else
    L = logL;
end

return;