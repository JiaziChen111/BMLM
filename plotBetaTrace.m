% Plot Beta 

betaAvg = squeeze(mean(beta_collection,1));

plot(betaAvg(1,:));
plot(betaAvg(2,:));
plot(betaAvg(3,:));
plot(betaAvg(4,:));