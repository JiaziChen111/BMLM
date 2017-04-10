function L = MyLogit(para, X1, X2, Q, flagLog)
% Likelihood of a multinomial logit function
% flagLog = 1, in log scale; 0, in normal scale
% Always set flagLog = 1 if calling MyLogit as the function to be optimized

% X_1:              J X D1 X T matrix
% X_2:              J X D2 matrix
% Q:                T X J matrix
% beta:             1 X D1 vector
% theta:            1 X D2 vector

T = size(X1, 3);
D1 = size(X1, 2);
D2 = size(X2, 2);

beta = para(1 : D1);
theta = para((D1 + 1) : (D1 + D2));

logL = 0;

if isempty(X2)
    for t = 1 : T
        a = exp( X1(:, :, t) * beta');
        logL = logL + sum( Q(t, : ) * ( log(a) - log(sum(a)) ) );
    end
else
    for t = 1 : T
        a = exp( X1(:, :, t) * beta' + X2 * theta');
        logL = logL + sum( Q(t, : ) * ( log(a) - log(sum(a)) ) );
    end
end

% scale
if flagLog == 1
    L = logL;
else
    L = exp( logL );
end

return;