function L = NormalLikelihood(x, mu, Sigma, invSigma, detSigma, flagFast, flagLog)
% Setting 1:
% Likelihood of X, a single sample, drawn from a normal with M sets of parameters
% x:                1 X D vector
% mu:               M X D vector
% Sigma:            D X D X M matrix
% invSigma:         D X D X M matrix (optional)
% detSigma:         M X 1 vector (optional)
% Setting 2:
% Likelihood of X, N samples, drawn from a normal. Each sample has a different correponding normal distribution
% x:                N X D vector
% mu:               N X D vector
% Sigma:            D X D X N matrix
% invSigma:         D X D X N matrix (optional)
% detSigma:         N X 1 vector (optional)
% Setting 3:
% Likelihood of X, N samples, drawn from a normal. All samples share the same normal distribution
% x:                N X D vector
% mu:               1 X D vector
% Sigma:            D X D matrix
% invSigma:         D X D matrix (optional)
% detSigma:         scalar (optional)
% flagFast: Scalar. 1 = use precomputed inverse and determination. invSigma
%                   and detSigma cannot be empty. 0 = use Sigma. (optional)
% flagLog:  Scalar. 1 = in log scale, 0 = in normal scale, default. (optional)
% Output:
% L:                M X 1 vector, or N X 1 vector

if nargin < 6
    flagFast = 0;
end

if nargin < 7
    flagLog = 0;
end

if flagFast == 0 && isempty(Sigma)
    fprintf('Error: Sigma cannot be empty!')
end

if flagFast == 1 && (isempty(invSigma) || isempty(detSigma))
    fprintf('Error: invSigma or detSigma cannot be empty!')
end

[N, D] = size(x);
M = size(mu, 1);

if N == 1 && M > 1 % if setting 1, format it the same as setting 2
    x = repmat(x, [M, 1]);
end

if N > 1 && M == 1 % if setting 3, format it the same as setting 2
    mu = repmat(mu, [N, 1]);
    Sigma = repmat(Sigma, [1, 1, N]);
    invSigma = repmat(invSigma, [1, 1, N]);
    detSigma = repmat(detSigma, [N, 1]);
    M = N;
end

logL = zeros(M, 1);

if flagFast == 0
    if isempty(detSigma)
        for k = 1 : M
            detSigma(k) = det(Sigma(:, :, k));
        end
    end
    df = x - mu; % M X D matrix
    logL = MultVMV(df, Sigma, df, 1); % df * inv(Sigma) * df'
    logL = - log(detSigma) / 2 - D / 2 * log(2 * pi) - logL / 2;
else
    df = x - mu; % M X D matrix
    logL = MultVMV(df, invSigma, df);
    logL = - log(detSigma) / 2 - D / 2 * log(2 * pi) - logL / 2;
end

if flagLog == 0
    L = exp( logL );
else
    L = logL;
end

return;