clear; clc; close all;
rng('default');

%% Test MultiVMV()
if 0
    D = 2;
    N = 10;
    v1 = rand(N, D);
    v2 = rand(N, D);
    Sigma = repmat(eye(D) * 5, [1,1,N]) + randn(D, D, N);
    % slow way
    for i = 1 : N
        xSlow(i, 1) = v1(i, :) * Sigma(:, :, i) * v2(i, :)';
    end
    xFast = MultVMV(v1, Sigma, v2);
    if N < 100
        [xSlow, xFast]
    end
    
    % slow way
    xSlow = [];
    for i = 1 : N
        xSlow(i, :) = v1(i, :) * Sigma(:, :, i);
    end
    xFast = MultVMV(v1, Sigma);
    if N < 100
        [xSlow, xFast]
    end
end

%% Specify data format
% X_1:              J X D1 X T X N matrix (in Matlab the 1st and 2nd subscripts refers rwo and column)
% X_2:              J X D2 matrix
% Q:                T X J X N matrix
% beta:     N X D1 matrix
% theta:    N X D2 matrix

%% Test NormalLikelihood()
if 0
    D = 2;
    M = 10000;
    x = rand(1, D);
    mu = rand(M, D);
    Sigma = repmat(eye(D) * 5, [1,1,M]) + randn(D, D, M);
    for k = 1 : M
        detSigma(k, 1) = det(Sigma(:, :, k));
        invSigma(:, : , k) = inv(Sigma(:, :, k));
    end
    LSlow = NormalLikelihood(x, mu, Sigma, invSigma, detSigma, 0, 1);
    LFast = NormalLikelihood(x, mu, Sigma, invSigma, detSigma, 1, 1);
    [LSlow, LFast]
    
    x = rand(M, D);
    LSlow = NormalLikelihood(x, mu, Sigma, invSigma, detSigma, 0, 1);
    LFast = NormalLikelihood(x, mu, Sigma, invSigma, detSigma, 1, 1);
    [LSlow, LFast]
    
    mu = mu(1, :);
    Sigma = Sigma(:, :, 1);
    invSigma = invSigma(:, :, 1);
    detSigma = detSigma(1);
    LSlow = NormalLikelihood(x, mu, Sigma, invSigma, detSigma, 0, 1);
    LFast = NormalLikelihood(x, mu, Sigma, invSigma, detSigma, 1, 1);
    [LSlow, LFast]
end

%% Test MyMvnrnd
if 0
    D = 2;
    N = 10000;
    b = zeros(N, D);
    S = [1 0.7; 0.7 1];
    Sigma = repmat(S, [1,1,N]);
    L = chol(S, 'lower');
    LRep = repmat(L, [1,1,N]);
    x = MyMvnrnd(b, Sigma, LRep);
    figure,
    scatter(x(:,1), x(:,2));
    axis([-4 4 -4 4]);
    
    x = MyMvnrnd(b, S, L);
    figure,
    scatter(x(:,1), x(:,2));
    axis([-4 4 -4 4]);
end

%% Test MNLLikelihood
if 0
    J = 6;
    T = 12;
    N = 10;
    D1 = 2;
    D2 = 5;
    
    X1 = rand(J, D1, T, N);
    X2 = rand(J, D2);
    Q = rand(T, J, N);
    beta = rand(N, D1);
    theta = rand(N, D2);
    
    MNLLikelihood(X1, X2, Q, beta, theta, 1, 1)
    MNLLikelihood(X1, X2, Q, beta, theta, 1, 0)
end

%% Test MyMvnrndPrec
if 1
    D = 2;
    N = 10000;
    b = zeros(N, D);
    % precision
    P = [1 0.7; 0.7 1];
    
    x = mvnrnd(b, inv(P));
    figure,
    scatter(x(:,1), x(:,2));
    axis([-4 4 -4 4]);
    mean(x)
    inv(cov(x))
    
    x = MyMvnrndPrec(b, P);
    figure,
    scatter(x(:,1), x(:,2));
    axis([-4 4 -4 4]);
    mean(x)
    inv(cov(x))
end


