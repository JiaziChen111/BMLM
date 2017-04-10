function [b, Sigma] = DrawNewComponent(b_0, upsilon_0, kappa_0, S_0, LS_0, Sigma_0, nd)

if nargin < 6
    nd = 1;
end

if nd > 1
    b_0 = repmat(b_0, [nd, 1]);
end

% Variance
if (0)
    if nargin < 5
        for i = 1 : nd
            Sigma(:, :, i) = iwishrnd(S_0, upsilon_0);
        end
    else % fast computation
        for i = 1 : nd
            Sigma(:, :, i) = iwishrnd(S_0, upsilon_0, LS_0);
        end
    end   
    % Mean
    b = mvnrnd(b_0, Sigma / kappa_0);
    
else
    % fix covariance
    Sigma = repmat(Sigma_0, [1, 1, nd]);
    % Mean
    b = mvnrnd(b_0, repmat(S_0, [1, 1, nd]));
end


