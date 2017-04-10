function k = DrawCFromOtherCs(N, compSize, ci, F)

% Remove c_i from the correpsonding component
% F is additional weight term

% Afterwards the compoent size should be 0 if c_i is a singleton
compSize(ci) = compSize(ci) - 1;

p = compSize / (N - 1);

if nargin == 4
    p = p .* F;
    p = p / sum(p); % p should sum to 1
end

% x = mnrnd(1, p);
% 
% k = find(x);


k = find(cumsum(p) > rand(1), 1, 'first');

return;