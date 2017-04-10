function d = GetDistance(X1, X2, m)
% Computer the L-m distance between X1 and X2. m = 2 only for now.
% Input:
% X1: N1 X p matrix
% X2: N2 X p matrix
% m:  Integer > 0
% Output
% d:  N1 X N2 matrix

[N1, p] = size(X1);
[N2, p2] = size(X2);

if p ~= p2
    fprintf('Error: X1 and X2 must have the same number of columns!\n')
    d = [];
    return;
end

if m == 2
    d = zeros(N1, N2);
    for i = 1 : N1
        d(i, :) = sqrt(sum((repmat(X1(i, :), N2, 1) - X2) .^ 2, 2))';
    end
end

return;

