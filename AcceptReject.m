function flag = AcceptReject(p)
% p:        Acceptance rate. Scalor ar N X 1 vector.

n = length(p);

x = rand(n, 1);

flag = x <= p;

return;

