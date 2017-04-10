function [b0, stdp] = FitLogit(X1, X2, Q)

% X_1:              J X D1 X T matrix
% X_2:              J X D2 matrix
% Q:                T X J matrix
% beta:             1 X D1 vector
% theta:            1 X D2 vector

% Determines the size of the X matrix
[J, D1, T] = size(X1);
D2 = size(X2, 2);
P = D1 + D2;

% X = [];
% Y = [];
% for t = 1 : T
%     X = [X; [X1(:, :, t), X2]];
%     Y = [Y; diag(Q(t, :))];
% end
% 
% b0 = mnrfit(X, Y);
% stdp = NaN;
% 
% return;

% Chooses random starting points
b00 = ones(1, D1 + D2) / 2;
b00(1) = -abs(b00(1));


% Finds the beta vector minimizing SSE
Tol = 1e-6;
OptimizationOptions = optimset('TolFun', Tol, 'TolX', Tol, 'Display', 'iter', 'MaxIter', 200000, 'MaxFunEvals', 500000, 'LargeScale', 'off', 'Algorithm', {'levenberg-marquardt',0.1});

tic
[b0,resnorm,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,Jacobian] = lsqnonlin(@(b) MyLogit(b, X1, X2, Q, 1), b00,[],[], OptimizationOptions);
toc

% Computes the SE of the parameters
Jacobian = full(Jacobian);
stdp = sqrt(diag(resnorm * inv(Jacobian' * Jacobian) / (T * J - P)));
