function L = NormalLikelihoodTheta(x, mu, Sigma, detSigma, flagLog)
% Likelihood of X, N samples, drawn from a normal. All samples share the same normal distribution
% Optimized for HIGH-dimension data

% x:                N X D vector
% mu:               1 X D vector
% Sigma:            D X D matrix
% detSigma:         scalar. Could be empty

% flagLog:          scalar. 1 = in log scale, 0 = in normal scale, default. 
% Output:
% L:                M X 1 vector, or N X 1 vector

%% Check input arguments
if nargin < 5
    fprintf('Error: All the input arguments have to be specified. detSigma could be empty though.')
end

if isempty(Sigma)
    fprintf('Error: Sigma cannot be empty if flagFast == 0 !')
end

[N, D] = size(x);
M = size(mu, 1);

if M ~= 1
    fprintf('Error: mu must be a row vector for setting 3!');
end

%% N samples from the same distribution
df = (x - repmat(mu, N, 1))'; % D X N
logL = Sigma \ df; % D X N
logL = sum(df .* logL, 1)'; % N X 1

if isempty(detSigma)
    logL = - log(det(Sigma)) / 2 - D / 2 * log(2 * pi) - logL / 2;
else
    logL = - log(detSigma) / 2 - D / 2 * log(2 * pi) - logL / 2;
end
    

%% Log or normal scale
if flagLog == 0
    L = exp( logL );
else
    L = logL;
end

return;