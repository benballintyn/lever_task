function [] = plotSession(mouse,day)
d = load(['processed_data/' mouse '_ReProcessedData.mat']); d=d.ProcessedData;
spstar = load(['analysis/spstar.mat']); spstar=spstar.spstar;
rstar = load(['analysis/rstar.mat']); rstar=rstar.rstar;

bestSP = spstar(d{day}.Reward_L,d{day}.TotalTrialsCompleted,d{day}.Reward_S,d{day}.NumPressRequired_S);
bestSP = ceil(bestSP);
bestR = rstar(d{day}.TotalTrialsCompleted,bestSP,d{day}.Reward_L,d{day}.Reward_S,d{day}.NumPressRequired_S);
aborted = find(isnan(d{day}.SideRewarded));
Srewarded = find(d{day}.SideRewarded == 0);
Lrewarded = find(d{day}.SideRewarded == 1);
totalLR = length(Lrewarded)*d{day}.Reward_L;
totalSR = length(Srewarded)*d{day}.Reward_S;
totalPress = sum(d{day}.LeverPressesByTrial);
Robs = (totalSR + totalLR)/totalPress;
baselineR = (d{day}.TotalTrialsCompleted*d{day}.Reward_S)/(d{day}.TotalTrialsCompleted*d{day}.NumPressRequired_S);

plot(d{day}.SideChosen); hold on;
plot(d{day}.LeverPressesByTrial);
scatter(aborted,d{day}.LeverPressesByTrial(aborted),'k*')
scatter(Srewarded,d{day}.LeverPressesByTrial(Srewarded),'b*')
scatter(Lrewarded,d{day}.LeverPressesByTrial(Lrewarded),'r*')
plot(1:d{day}.TotalTrialsCompleted,ones(1,d{day}.TotalTrialsCompleted)*bestSP,'linew',2)
title(['R^* = ' num2str(bestR) ' R^{obs} = ' num2str(Robs) ' R_0 = ' num2str(baselineR)])
end

