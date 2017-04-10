function [X1, X2, Q] = LoadData(fileName, dimPara)
% Load data for MNL-DPM

% Unpack parameters
[N, T, J, D1, D2, NT, TJ, NTJ] = UnpackDimPara(dimPara);

% Load data
data = load(fileName);
if N ~=  size(data,1) / TJ
    fprintf('Error: input data size is not correct!');
end

% Q: #visits to a store
for i = 1 : N
    Q(:, :, i) = reshape(data((i - 1) * TJ + (1 : TJ), 5), J, T)';
end

% X1: 
X1 = zeros(J, D1, T, N);
for i = 1 : N
    for t = 1 : T
        X1(:, 1, t, i) = prod(data((i - 1) * TJ + (t - 1) * J + (1 : J), 7 : 8), 2); % price*distance
        if D1 > 1
            X1(:, 2, t, i) = data((i - 1) * TJ + (t - 1) * J + (1 : J), 8); %distance
        end
        if D1 > 2
            X1(:, 3, t, i) = data((i - 1) * TJ + (t - 1) * J + (1 : J), 7); %price
        end
    end
end

% X2
X2 = eye(J - 1); % identity matrix of dim J
X2 = [X2; zeros(1, J - 1)]; % "Other" stores

% Normalize X1
X1New = zeros(J, D1, T, N);
for i = 1 : N
    for t = 1 : T
        X1New(:, :,  t, i) = X1(:, :, t, i) - repmat(X1(J, :, t, i), J, 1);
    end
end
X1 = X1New;

return;

