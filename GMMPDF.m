function GMMPDF(w, b, Sigma)

x = -20 : 0.1 : 20;

[K, D] = size(b);

for j = 1 : D
    p = zeros(1, length(x));
    for k = 1 : K
        p = p + w(k) * pdf('norm', x, b(k, j), sqrt(Sigma(j, j, k)));
    end
     p = p / sum(p);
     
     figure,
     plot(x, p);
     title(['var ', num2str(j)]);
end