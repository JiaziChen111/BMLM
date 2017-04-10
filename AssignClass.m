function c = AssignClass(X, classCenter) 

D = GetDistance(X, classCenter, 2);

[~, c] = min(D, [], 2);

return;



