function [m, prRL, prG] = EvalSamp(samples, thinFactor, ids, ide, flagCoda, flagKDE, vName)

% here we assume the last subscript is for sample indices. e.g.
% beta(2,3,1000) has 1000 samples

s = size(samples);
% single sample
if s(end) == 1
    return;
end
d = length(s);
nSamp = s(end);

% thinning
if isempty(ids)
    ids = 1;
end
if isempty(ide)
    ide = nSamp;
end
idx = ids : thinFactor : ide;
switch d
    case 2
        samples = samples(:, idx);
    case 3
        samples = samples(:, :, idx);
    otherwise
        return;
end
nSamp = length(idx);

% plot KDE (over the first dimension) of the last sample
if flagKDE && s(1) > 1
    switch d
        case 2
            [f, xi] = ksdensity(mean(samples, 2), ...
                'function', 'pdf', 'width', 1);
            figure,
            plot(xi, f);
            title('kde for var 1');
        case 3
            for j = 1 : s(2)
                [f, xi] = ksdensity(mean(samples(:,j,:), 3), ...
                    'function', 'pdf', 'width', 1);
                figure,
                plot(xi, f);
                title(['kde for var ', num2str(j)]);
            end
            if s(2) == 2 %bivariate
                [bandwidth,density,X,Y] = kde2d(mean(samples, 3));
                figure,
                surf(X,Y,density);
                %hold on
                %plot(data(:,1),data(:,2),'r.','MarkerSize',5);
                axis([-6 6 -4 4]);
            end
        otherwise
            return;
    end
end

% sample mean
switch d
    case 2
        m = mean(samples, 2);
    case 3
        m = mean(samples, 3);
    otherwise
        return;
end

% plot sample chain for C
if flagKDE == 0 && length(s) == 2
    figure(101),
    plot(samples);
    title('chain');
end

% coda
if flagCoda
    switch d
        case 2
            x = samples';
        case 3
            x = [];
            for i = 1 : s(1)
                beta = squeeze(samples(i,:,:))';
                if size(beta, 1) == 1
                    beta = beta';
                end
                x = [x, beta];
            end
        otherwise
            return;
    end
    r = coda(x);
    nVar = size(x, 2);
    for j = 1 : nVar
        IScore(j) = r(j).irl;
        pchisqr(j) = max(r(j).pchisqr);
    end
    % convergency-check pass rate
    prRL = sum(IScore < 5) / nVar;
    prG = sum(pchisqr > 0.05) / nVar;
    
    % print out
    fprintf([vName, ':\n']);
    %fprintf('mean = %f\n', m);
    fprintf('convergency RL = %f\n', prRL);
    fprintf('convergency G = %f\n', prG);
    
else
    prRL = [];
    prG = [];
end