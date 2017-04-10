function L = NormalLikelihoodBeta(x, mu, Sigma, invSigma, detSigma, setting, flagFast, flagLog)
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
% invSigma:         D X D matrix 
% detSigma:         scalar 

% setting:          scalar, 1, 2 or 3.
% flagFast:         scalar. 1 = use precomputed inverse and determination. invSigma
%                   and detSigma cannot be empty. 0 = use Sigma. (optional)
% flagLog:          scalar. 1 = in log scale, 0 = in normal scale, default. 
% Output:
% L:                M X 1 vector, or N X 1 vector

%% Check input arguments
if nargin < 8
    fprintf('Error: All the input arguments have to be specified. Sigma/invSigma/detSigma could be empty though.')
end

if flagFast == 0 && isempty(Sigma)
    fprintf('Error: Sigma cannot be empty if flagFast == 0 !')
end

if flagFast == 1 && (isempty(invSigma) || isempty(detSigma))
    fprintf('Error: invSigma or detSigma cannot be empty if flagFast == 1!')
end

[N, D] = size(x);
M = size(mu, 1);

if setting == 1 && N ~= 1
    fprintf('Error: x must be a row vector for setting 1!');
end

if setting == 2 && M ~= N
    fprintf('Error: x must have the same size as mu for setting 2!');
end

if setting == 3 && M ~= 1
    fprintf('Error: mu must be a row vector for setting 3!');
end

%% N samples from the same distribution
if setting == 3 
    if flagFast
        df = (x - repmat(mu, N, 1))'; % D X N
        logL = invSigma * df; % D X N
        logL = sum(df .* logL, 1)'; % N X 1
        logL = - log(detSigma) / 2 - D / 2 * log(2 * pi) - logL / 2;
    else
        df = (x - repmat(mu, N, 1))'; % D X N
        logL = Sigma \ df; % D X N
        logL = sum(df .* logL, 1)'; % N X 1
        logL = - log(det(Sigma)) / 2 - D / 2 * log(2 * pi) - logL / 2;
    end
end

%% M samples from M different distributions
if setting == 1 % if setting 1, format it the same as setting 2
    x = repmat(x, [M, 1]);
    N = M;
end

if isempty(detSigma)
    for k = 1 : M
        detSigma(k) = det(Sigma(:, :, k));
    end
end

df = x - mu; % M X D matrix
logL = zeros(M, 1);

% batches for big M(N)
batch = 1000;
nBatch = ceil(N / batch);
for iB = 1 : nBatch
    if iB < nBatch
        i = (iB - 1) * batch + 1 : iB * batch;
    else
        i = (iB - 1) * batch + 1 : N;
    end
    if flagFast
        logL(i) = MultVMV(df(i, :), invSigma(:, :, i), df(i, :));
    else
        logL(i) = MultVMV(df(i, :), Sigma(:, :, i), df(i, :), 1); % df * inv(Sigma) * df'
    end
end

logL = - log(detSigma) / 2 - D / 2 * log(2 * pi) - logL / 2;

%% Log or normal scale
if flagLog == 0
    L = exp( logL );
else
    L = logL;
end

return;