function x = MyMvnrnd(b, Sigma, L)
% Setting 1:
% b:                1 X D
% Sigma:            D X D
% L:                D X D
% Setting 2:
% b:                N X D 
% Sigma:            D X D
% L:                D X D
% Setting 3:
% b:                N X D 
% Sigma:            D X D X N
% L:                D X D X N


if nargin < 3
    x = mvnrnd(b, Sigma);
else
    [N, D] = size(b);
    if length(size(L)) == 2 && N > 1 %setting 2
        LT = repmat(L', [1, 1, N]);
    else 
        LT = permute(L, [2, 1, 3]);
    end
    x = randn(N, D);
    x = MultVMV(x, LT); % 
    x = x + b;
%     x = repmat(reshape(randn(D*N, 1), [1, D, N]), [D, 1, 1]); % D X D X N
%     if D == 1
%         x = squeeze(sum(L .* x, 2)) + b; % N X 1
%     else
%         x = squeeze(sum(L .* x, 2))' + b; % N X D
%     end
end

return;