% make pizza talk figures
NL_optimal = [10 22 28 58];
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
NL_observed = load(['optimality/WT/NL_observed.mat']); NL_observed=NL_observed.NL_observed;
perfectAgentNumLR = load('pizza_talk_2020/data_files/perfectAgentNumLR.mat'); perfectAgentNumLR=perfectAgentNumLR.perfectAgentNumLR;
figFolder = 'pizza_talk_2020/figures/';

%% Figure 1 Histogram of NL_observed with NL_optimal and perfectAgentNumLR labled
figure;
for i = 1:4
    lineHeight = sum(NL_observed{i} == mode(NL_observed{i}));
    subplot(2,2,i);
    histogram(NL_observed{i},'binWidth',1)
    hold on;
    plot([NL_optimal(i) NL_optimal(i)],[0 lineHeight],'Color','k')
    xlabel('# of PR trials')
    ylabel('# of sessions')
    if (i == 1)
        legend({'Mice data','Optimal'})
    end
    xlim([0 90])
    set(gca,'fontsize',15,'fontweight','bold')
end
set(gcf,'Position',[10 10 1600 1200])
saveas(gcf,[figFolder 'NL_observed_NL_optimal_perfectAgentNumLR.fig'],'fig')
saveas(gcf,[figFolder 'NL_observed_NL_optimal_perfectAgentNumLR.eps'],'eps')
saveas(gcf,[figFolder 'NL_observed_NL_optimal_perfectAgentNumLR.tif'],'tiffn')

%% Perfect Agent Q-values vs. trial number with optimal label
load('~/phd/lever_task/pizza_talk_2020/data_files/perfectAgentQ.mat');
for i=1:4
    subplot(2,2,i)
    plot(Q{i})
    hold on;
    plot([NL_optimal(i) NL_optimal(i)],[0 8],'Color','k')
    legend({'FR','PR','Optimal PR #'})
    title(sessionTypes{i})
    ylim([0 8])
    xlabel('PR trial #','fontsize',20,'fontweight','bold')
    ylabel('EoR','fontsize',20,'fontweight','bold')
    set(gca,'fontsize',15,'fontweight','bold')
end
set(gcf,'Position',[10 10 1600 1200])
saveas(gcf,[figFolder 'perfect_agent_Qvalues_optimal_labeled.fig'],'fig')
saveas(gcf,[figFolder 'perfect_agent_Qvalues_optimal_labeled.eps'],'eps')
saveas(gcf,[figFolder 'perfect_agent_Qvalues_optimal_labeled.tif'],'tiffn')

%% Random corrected EoRs (Paul method)
EoR_optimalities = load(['~/phd/lever_task/optimality/WT/EoR_optimalities.mat']); EoR_optimalities = EoR_optimalities.EoR_optimalities;
EoRRandOptimalities = load(['~/phd/lever_task/optimality/WT/EoRRandOptimalities.mat']); EoRRandOptimalities = EoRRandOptimalities.EoRRandOptimalities;
nTrials = load(['~/phd/lever_task/optimality/WT/nTrials.mat']); nTrials=nTrials.nTrials;
for i=1:4
    for j=1:length(EoR_optimalities{i})
        correctedEoRs{i}(j) = (EoR_optimalities{i}(j) - mean(EoRRandOptimalities{i}(nTrials{i}(j),:)))/(1 - mean(EoRRandOptimalities{i}(nTrials{i}(j),:)));
    end
    correctedEoRsNoOutliers{i} = correctedEoRs{i}(correctedEoRs{i} > -5);
    disp([sessionTypes{i} ': mean +- SD = ' num2str(mean(correctedEoRs{i})) ' +- ' num2str(std(correctedEoRs{i}))]) 
end
save('~/phd/lever_task/pizza_talk_2020/data_files/correctedEoRs.mat','correctedEoRs','-mat')
save('~/phd/lever_task/optimality/WT/correctedEoRs.mat','correctedEoRs','-mat')

figure;
for i=1:4
    notBoxPlot(correctedEoRs{i},i)
end
set(gca,'xtick',1:4,'xticklabels',{'2xFR6','2xFR12','5xFR6','5xFR12'})
xlabel('Session Type'); ylabel('Corrected EoR optimality')
saveas(gcf,[figFolder 'correctedEoRs.fig'],'fig')
saveas(gcf,[figFolder 'correctedEoRs.eps'],'eps')
saveas(gcf,[figFolder 'correctedEoRs.tif'],'tif')

figure;
for i=1:3
    notBoxPlot(correctedEoRs{i},i)
end
set(gca,'xtick',1:3,'xticklabels',{'2xFR6','2xFR12','5xFR6'})
xlabel('Session Type'); ylabel('Corrected EoR optimality')
saveas(gcf,[figFolder 'correctedEoRs_2xFR6_2xFR12_5xFR6_only.fig'],'fig')
saveas(gcf,[figFolder 'correctedEoRs_2xFR6_2xFR12_5xFR6_only.eps'],'eps')
saveas(gcf,[figFolder 'correctedEoRs_2xFR6_2xFR12_5xFR6_only.tif'],'tif')

figure;
for i=1:4
    notBoxPlot(correctedEoRs{i},i)
end
set(gca,'xtick',1:4,'xticklabels',{'2xFR6','2xFR12','5xFR6','5xFR12'})
ylim([-3 1])
xlabel('Session Type'); ylabel('Corrected EoR optimality')
saveas(gcf,[figFolder 'correctedEoRs_zoom.fig'],'fig')
saveas(gcf,[figFolder 'correctedEoRs_zoom.eps'],'eps')
saveas(gcf,[figFolder 'correctedEoRs_zoom.tif'],'tif')

figure;
for i=1:4
    notBoxPlot(correctedEoRsNoOutliers{i},i);
end
set(gca,'xtick',1:4,'xticklabels',{'2xFR6','2xFR12','5xFR6','5xFR12'})
xlabel('Session Type'); ylabel('Corrected EoR optimality')
saveas(gcf,[figFolder 'correctedEoRs_noOutliers.fig'],'fig')
saveas(gcf,[figFolder 'correctedEoRs_noOutliers.eps'],'eps')
saveas(gcf,[figFolder 'correctedEoRs_noOutliers.tif'],'tif')

correctedEoRStats_WT

%% Show computation of corrected EoR
EoR_optimalities = load(['~/phd/lever_task/optimality/WT/EoR_optimalities.mat']); EoR_optimalities = EoR_optimalities.EoR_optimalities;
EoRRandOptimalities = load(['~/phd/lever_task/optimality/WT/EoRRandOptimalities.mat']); EoRRandOptimalities = EoRRandOptimalities.EoRRandOptimalities;
nTrials = load(['~/phd/lever_task/optimality/WT/nTrials.mat']); nTrials=nTrials.nTrials;

max5xFR6_ntrial = find(max(nTrials{3}));
figure;
plot(mean(EoRRandOptimalities{3},2));
hold on;
scatter(nTrials{3}(max5xFR6_ntrial),EoR_optimalities{3}(max5xFR6_ntrial),500,'.')
xlabel('Trial #')
ylabel('EoR optimality')
set(gca,'fontsize',15,'fontweight','bold')
legend({'EoR optimality of random agent','EoR optimality of mouse'},'location','southwest')
saveas(gcf,[figFolder 'corrected_EoR_computation_example.fig'],'fig')
saveas(gcf,[figFolder 'corrected_EoR_computation_example.eps'],'eps')
saveas(gcf,[figFolder 'corrected_EoR_computation_example.tif'],'tiffn')

%% Zoomed in plot of P(FR) and P(PR)
pS = load('~/phd/lever_task/optimality/WT/pS.mat'); pS=pS.pS;
pL = load('~/phd/lever_task/optimality/WT/pL.mat'); pL=pL.pL;
figure;
for i=1:4
    subplot(2,2,i)
    plot(pS{i}(1:300)); hold on; plot(pL{i}(1:300))
    xlabel('Trial number')
    ylabel('Probability')
    legend({'P(FR)','P(PR)'})
    title(sessionTypes{i})
    set(gca,'fontsize',15,'fontweight','bold')
end
set(gcf,'Position',[10 10 1600 1200])
saveas(gcf,[figFolder 'p_FR_p_PR_zoom.fig'],'fig')
saveas(gcf,[figFolder 'p_FR_p_PR_zoom.eps'],'eps')
saveas(gcf,[figFolder 'p_FR_p_PR_zoom.tif'],'tiffn')

%% Plot of P(FR) and P(PR) from best match agent
bestMatchInds = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_ansUtilityFunc_pressTimeFunc_initialization_meanReward_ans0-1_analyzed/bestMatchInds.mat');
bestMatchInds = bestMatchInds.bestMatchInds;

pS = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_ansUtilityFunc_pressTimeFunc_initialization_meanReward_ans0-1_analyzed/pS.mat');
pS = pS.pS;

pL = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_ansUtilityFunc_pressTimeFunc_initialization_meanReward_ans0-1_analyzed/pL.mat');
pL = pL.pL;

figure;
for i=1:4
    subplot(2,2,i)
    plot(pS{bestMatchInds(i,1),bestMatchInds(i,2),bestMatchInds(i,3),i}(1:300));
    hold on;
    plot(pL{bestMatchInds(i,1),bestMatchInds(i,2),bestMatchInds(i,3),i}(1:300));
    if (i == 1)
        legend({'P(FR)','P(PR)'})
    end
    if (i == 1 || i == 3)
        ylabel('Probability')
    end
    if (i == 3 || i == 4)
        xlabel('Trial #')
    end
    set(gca,'fontsize',15,'fontweight','bold')
end
set(gcf,'Position',[10 10 1600 1200])
saveas(gcf,[figFolder 'agent_p_FR_p_PR_zoom.fig'],'fig')
saveas(gcf,[figFolder 'agent_p_FR_p_PR_zoom.eps'],'eps')
saveas(gcf,[figFolder 'agent_p_FR_p_PR_zoom.tif'],'tiffn')

%% Give example of approximate number system
sigmas = [.1 .25 .4 .7];
figure;
for i=1:4
    subplot(2,2,i)
    histogram(lognrnd(log(50) - (sigmas(i)^2)/2,sigmas(i),1,4000),'normalization','pdf','binwidth',1)
    if (i == 3 || i == 4)
        xlabel('Estimation of number 50')
    end
    if (i == 1 || i == 3)
        ylabel('Probability density')
    end
    hold on;
    lineHeight = ylim;
    plot([50 50],[0 lineHeight(2)])
end
set(gcf,'Position',[10 10 1400 1200])
saveas(gcf,[figFolder 'ANS_example.fig'],'fig')
saveas(gcf,[figFolder 'ANS_example.eps'],'eps')
saveas(gcf,[figFolder 'ANS_example.tif'],'tiffn')

%% Show how increasing ANS sigma shifts PR trial distribution
bestNumPR = load(['~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_ansUtilityFunc_pressTimeFunc_initialization_meanReward_ans0-1_analyzed/numLR.mat']);
bestNumPR = bestNumPR.numLR;
bestMatchInds = load(['~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_ansUtilityFunc_pressTimeFunc_initialization_meanReward_ans0-1_analyzed/bestMatchInds.mat']);
bestMatchInds = bestMatchInds.bestMatchInds;
figure;
histogram(bestNumPR{bestMatchInds(4,1),bestMatchInds(4,2),5,4},'normalization','pdf','binwidth',1)
hold on;
histogram(bestNumPR{bestMatchInds(4,1),bestMatchInds(4,2),10,4},'normalization','pdf','binwidth',1)
%histogram(bestNumPR{bestMatchInds(4,1),bestMatchInds(4,2),20,4},'normalization','pdf','binwidth',1)
xlabel('PR trials completed','fontsize',20,'fontweight','bold')
ylabel('Probability Density','fontsize',20,'fontweight','bold')
legend({'ANS \sigma = .25','ANS \sigma = .5'})
set(gca,'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1000])
saveas(gcf,[figFolder 'PR_trials_completed_ANS_sigma_shift.fig'],'fig')
saveas(gcf,[figFolder 'PR_trials_completed_ANS_sigma_shift.eps'],'eps')
saveas(gcf,[figFolder 'PR_trials_completed_ANS_sigma_shift.tif'],'tiffn')

%% Softmax example
vals = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1 .9 .8 .7 .6 .5 .4 .3 .2 .1];
[softvals,~] = mySoftmax(vals,1);
[softvals2,~] = mySoftmax(vals,.5);
[softvals3,~] = mySoftmax(vals,.2);
figure;
scatter(vals,softvals,50,'filled'); hold on; 
scatter(vals,softvals2,50,'filled')
scatter(vals,softvals3,50,'filled')
legend({'\tau = 1','\tau = .5','\tau = .2'})
xlabel('Q-value of action')
ylabel('Probability of action')
saveas(gcf,[figFolder 'softmax_example.fig'],'fig')
saveas(gcf,[figFolder 'softmax_example.eps'],'eps')
saveas(gcf,[figFolder 'softmax_example.tif'],'tiffn')

%% Softmax as logistic function example
x = -.5:.01:5;
f = @(x,tau) 1./(1 + exp(-(1/tau)*x));
plot(x,f(x,.2)); hold on; plot(x,f(x,.05));
set(gca,'fontsize',15,'fontweight','bold')
legend({'\tau = .1','\tau = .05'},'fontsize',30,'fontweight','bold','location','east');
xlabel('Q(PR) - Q(FR)')
ylabel('P(PR)');
saveas(gcf,[figFolder 'softmax_as_logistic_fxn_example.fig'],'fig')
saveas(gcf,[figFolder 'softmax_as_logistic_fxn_example.eps'],'eps')
saveas(gcf,[figFolder 'softmax_as_logistic_fxn_example.tif'],'tiffn')

%% Comparison of fits across agent types
wd_softmax = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_initialization_trained_analyzed/wasserstein_distances.mat');
wd_softmax = wd_softmax.wasserstein_distances;
bestMatchInds_softmax = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_initialization_trained_analyzed/bestMatchInds.mat');
bestMatchInds_softmax = bestMatchInds_softmax.bestMatchInds;

wd_softmax_pressTimeFunc = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_pressTimeFunc_initialization_trained_analyzed/wasserstein_distances.mat');
wd_softmax_pressTimeFunc = wd_softmax_pressTimeFunc.wasserstein_distances;
bestMatchInds_softmax_pressTimeFunc = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_pressTimeFunc_initialization_trained_analyzed/bestMatchInds.mat');
bestMatchInds_softmax_pressTimeFunc = bestMatchInds_softmax_pressTimeFunc.bestMatchInds;

wd_softmax_ansUtilityFunc = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_ansUtilityFunc_initialization_trained_analyzed/wasserstein_distances.mat');
wd_softmax_ansUtilityFunc = wd_softmax_ansUtilityFunc.wasserstein_distances;
bestMatchInds_softmax_ansUtilityFunc = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_ansUtilityFunc_initialization_trained_analyzed/bestMatchInds.mat');
bestMatchInds_softmax_ansUtilityFunc = bestMatchInds_softmax_ansUtilityFunc.bestMatchInds;

wd_softmax_ansUtilityFunc_pressTimeFunc = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_ansUtilityFunc_pressTimeFunc_initialization_trained_analyzed/wasserstein_distances.mat');
wd_softmax_ansUtilityFunc_pressTimeFunc = wd_softmax_ansUtilityFunc_pressTimeFunc.wasserstein_distances;
bestMatchInds_softmax_ansUtilityFunc_pressTimeFunc = load('~/phd/lever_task/simple_rl/results/cluster_results/bandit_softmax_ansUtilityFunc_pressTimeFunc_initialization_trained_analyzed/bestMatchInds.mat');
bestMatchInds_softmax_ansUtilityFunc_pressTimeFunc = bestMatchInds_softmax_ansUtilityFunc_pressTimeFunc.bestMatchInds;
figure;
for i=1:4
    subplot(2,2,i)
    notBoxPlot(wd_softmax(bestMatchInds_softmax(i,1),bestMatchInds_softmax(i,2)),1)
    hold on;
    notBoxPlot(wd_softmax_pressTimeFunc(bestMatchInds_softmax_pressTimeFunc(i,1),bestMatchInds_softmax_pressTimeFunc(i,2)),2)
    notBoxPlot(wd_softmax_ansUtilityFunc(bestMatchInds_softmax_ansUtilityFunc(i,1),bestMatchInds_softmax_ansUtilityFunc(i,2),bestMatchInds_softmax_ansUtilityFunc(i,3)),3)
    notBoxPlot(wd_softmax_ansUtilityFunc_pressTimeFunc(bestMatchInds_softmax_ansUtilityFunc_pressTimeFunc(i,1),bestMatchInds_softmax_ansUtilityFunc_pressTimeFunc(i,2),bestMatchInds_softmax_ansUtilityFunc_pressTimeFunc(i,3)),4)
end