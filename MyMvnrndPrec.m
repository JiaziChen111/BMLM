function x = MyMvnrndPrec(B, Lambda)
% Function to draw n random samples from a normal distribution given mean b
% and precision Lambda
% Use this function when inverse of the precision is difficult, e.g. large matrix
% B:          N X D vector
% Lambda:     D X D matrix
% n:          scalar

[N, D] = size(B);

R = chol(Lambda); % upper triangle matrix

v = randn(D, N);
x = (R\v)' + B;

% Proof: Covariance of X (n X D matrix):
% (X - B)' * (X - B) = (R\V) * (R\V)' = inv(R) * V * V' * inv(R)'
% = inv(R) * inv(R)' = inv(R' * R) = inv(Lambda)

return;