function [prod1, prod2] = PrecProd(X_1, X_2, beta, theta, N, T, J)

% X_1:              (N*TJ) X (N*D1) sparse matrix 
% X_2:              J X D2 matrix
% beta:             N X D1 matrix
% theta:            N X D2 matrix

if N > 1
    % prod1 = squeeze(sum(X_1 .* repmat(reshape(beta', [1, D1, 1, N]), [J, 1, T, 1]), 2)); % J X T X N
    prod1 = reshape(X_1 * reshape(beta', [], 1), [J, T, N]);
    prod2 = repmat(reshape([theta'; zeros(1, N)], [J, 1, N]), [1, T, 1]); % J X T X N
else
    % prod1 = squeeze(sum(X_1 .* repmat(beta, [J, 1, T]), 2));
    prod1 = reshape(X_1 * beta, [J, T]);
    prod2 = repmat([theta'; 0], [1, T]);
end

return;