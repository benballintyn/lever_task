%pizza_talk_figs_2021
%%
figureDir = '/home/ben/phd/lever_task/pizza_talk_2021/';
if (~exist(figureDir,'dir'))
    mkdir(figureDir)
end

useOnly120Trials = 1;

externalHDDir = '/media/ben/Varda/';
externalDataDir = [externalHDDir 'phd/lever_task/'];
if (useOnly120Trials)
    driftRLDir = [externalDataDir 'driftRL/results/Only120Trials/logprob_joint/value_based_drift/'];
else
    driftRLDir = [externalDataDir 'driftRL/results/logprob_joint/value_based_drift/'];
end
driftSimDataDirs{1} = 'bandit_e_greedy_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
driftSimDataDirs{2} = 'bandit_e_greedy_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
driftSimDataDirs{3} = 'bandit_e_greedy_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
driftSimDataDirs{4} = 'bandit_e_greedy_initialization_mean_reward_forgettingType_decayToInitialValues';
driftSimDataDirs{5} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
for i=1:length(driftSimDataDirs)
    driftRLDirs{i} = [driftRLDir driftSimDataDirs{i}];
end

if (useOnly120Trials)
    logisticRLDir = [externalDataDir 'logisticAbortRL/results/Only120Trials/logprob_joint/'];
else
    logisticRLDir = [externalDataDir 'logisticAbortRL/results/logprob_joint/'];
end
logisticSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
logisticSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
logisticSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
logisticSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues';
logisticSimDataDirs{5} = 'bandit_e_greedy_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';

for i=1:length(logisticSimDataDirs)
    logisticRLDirs{i} = [logisticRLDir logisticSimDataDirs{i}];
end

if (useOnly120Trials)
    driftRL_valueUpdateDir = [externalDataDir 'driftRL_valueUpdate/results/Only120Trials/logprob_joint/value_based_drift/'];
else
    driftRL_valueUpdateDir = [externalDataDir 'driftRL_valueUpdate/results/logprob_joint/value_based_drift/'];
end
driftRL_valueUpdateSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
driftRL_valueUpdateSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
driftRL_valueUpdateSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
driftRL_valueUpdateSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues';
for i=1:length(driftRL_valueUpdateSimDataDirs)
    driftRL_valueUpdateDirs{i} = [driftRL_valueUpdateDir driftRL_valueUpdateSimDataDirs{i}];
end

%% DriftRL figs first
for i=1:length(driftRLDirs)
    disp(['Making plots for: ' driftRLDirs{i}])
    [scores,allParams] = plotResults(driftRLDirs{i},1,'fitStatsOnly',true);
    saveas(gcf,[figureDir 'driftRL_' driftSimDataDirs{i} '_bestFitStats.fig'],'fig')
    print([figureDir 'driftRL_' driftSimDataDirs{i} '_bestFitStats.png'],'-dpng','-r600')
end

%% LogisticAbort figs
for i=1:length(logisticRLDirs)
    disp(['Making plots for: ' logisticRLDirs{i}])
    [scores,allParams] = plotResults(logisticRLDirs{i},1,'fitStatsOnly',true);
    saveas(gcf,[figureDir 'logisticAbortRL_' logisticSimDataDirs{i} '_bestFitStats.fig'],'fig')
    print([figureDir 'logisticAbortRL_' logisticSimDataDirs{i} '_bestFitStats.png'],'-dpng','-r600')
end

%% driftRL_valueUpdate figs
for i=1:length(driftRL_valueUpdateDirs)
    disp(['Making plots for: ' driftRL_valueUpdateDirs{i}])
    [scores,allParams] = plotResults(driftRL_valueUpdateDirs{i},1,'fitStatsOnly',true);
    saveas(gcf,[figureDir 'driftRL_valueUpdate_' driftRL_valueUpdateSimDataDirs{i} '_bestFitStats.fig'],'fig');
    print([figureDir 'driftRL_valueUpdate_' driftRL_valueUpdateSimDataDirs{i} '_bestFitStats.png'],'-dpng','-r600')
end

%% Aborted trials analysis
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
nTrialsAborted = load('optimality/WT/nTrialsAborted.mat'); nTrialsAborted=nTrialsAborted.nTrialsAborted;
nTrialsCompleted = load('optimality/WT/nTrialsCompleted.mat'); nTrialsCompleted=nTrialsCompleted.nTrialsCompleted
fracTrialsAborted = cell(1,4);
for i=1:4
    for j=1:length(nTrialsCompleted{i})
        fracTrialsAborted{i} = [fracTrialsAborted{i} nTrialsAborted{i}(j)/(nTrialsAborted{i}(j) + nTrialsCompleted{i}(j))];
    end
end
figure;
for i=1:4
    subplot(2,2,i)
    histogram(fracTrialsAborted{i},'normalization','pdf','binwidth',.02)
    title(sessionTypes{i},'FontSize',15,'FontWeight','bold')
    xlabel('Fraction of trials aborted','FontSize',15,'FontWeight','bold')
    if (i == 1 || i == 3)
        ylabel('Probability density','FontSize',15,'FontWeight','bold')
    end
    xlim([0 .5])
end
set(gcf,'Position',[10 10 1600 1600])
saveas(gcf,[figureDir 'fracTrialsAborted.fig'],'fig')
print([figureDir 'fracTrialsAborted.png'],'-dpng','-r600')

%% % PR trial completed closeup
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
percentCompletedPR_allTrials = load('driftRL/abort_analysis/percentCompletedPR_allTrials.mat'); percentCompletedPR_allTrials = percentCompletedPR_allTrials.percentCompletedPR_allTrials;
% Logistic RL first
[scores,allParams] = plotlogisticAbortRLResults(logisticRLDirs{1},1,'fitStatsOnly',true);
[~,maxInd] = min(scores);
percentCompletedPR = load([logisticRLDirs{1} '/' num2str(maxInd) '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
figure;
for i=1:4
    subplot(2,2,i)
    curvals = percentCompletedPR_allTrials{i}(percentCompletedPR_allTrials{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    hold on;
    curvals = percentCompletedPR{i}(percentCompletedPR{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    title(sessionTypes{i},'fontsize',15,'fontweight','bold')
    xlabel('% of PR trial completed','fontsize',15,'fontweight','bold')
    if (i == 1)
        legend({'Mice','Model'},'fontsize',15,'fontweight','bold','Location','northeast')
    end
    if (i == 1 || i == 3)
        ylabel('Probability density','fontsize',15,'fontweight','bold')
    end
end
set(gcf,'Position',[10 10 1600 1600])
saveas(gcf,[figureDir 'percentCompletedPR_zoom_logisticRL.fig'],'fig')
print([figureDir 'percentCompletedPR_zoom_logisticRL.png'],'-dpng','-r600')

% Drift RL next
[scores,allParams] = plotlogisticAbortRLResults(driftRLDirs{1},1,'fitStatsOnly',true);
[~,maxInd] = min(scores);
percentCompletedPR = load([driftRLDirs{1} '/' num2str(maxInd) '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
figure;
for i=1:4
    subplot(2,2,i)
    curvals = percentCompletedPR_allTrials{i}(percentCompletedPR_allTrials{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    hold on;
    curvals = percentCompletedPR{i}(percentCompletedPR{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    title(sessionTypes{i},'fontsize',15,'fontweight','bold')
    xlabel('% of PR trial completed','fontsize',15,'fontweight','bold')
    if (i == 1)
        legend({'Mice','Model'},'fontsize',15,'fontweight','bold','Location','northeast')
    end
    if (i == 1 || i == 3)
        ylabel('Probability density','fontsize',15,'fontweight','bold')
    end
end
set(gcf,'Position',[10 10 1600 1600])
saveas(gcf,[figureDir 'percentCompletedPR_zoom_driftRL.fig'],'fig')
print([figureDir 'percentCompletedPR_zoom_driftRL.png'],'-dpng','-r600')

% DriftRL_valueUpdate last
[scores,allParams] = plotDriftRL_valueUpdateResults(driftRL_valueUpdateDirs{1},1,'fitStatsOnly',true);
[~,maxInd] = min(scores);
percentCompletedPR = load([driftRL_valueUpdateDirs{1} '/' num2str(maxInd) '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
figure;
for i=1:4
    subplot(2,2,i)
    curvals = percentCompletedPR_allTrials{i}(percentCompletedPR_allTrials{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    hold on;
    curvals = percentCompletedPR{i}(percentCompletedPR{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    title(sessionTypes{i},'fontsize',15,'fontweight','bold')
    xlabel('% of PR trial completed','fontsize',15,'fontweight','bold')
    if (i == 1)
        legend({'Mice','Model'},'fontsize',15,'fontweight','bold','Location','northeast')
    end
    if (i == 1 || i == 3)
        ylabel('Probability density','fontsize',15,'fontweight','bold')
    end
end
set(gcf,'Position',[10 10 1600 1600])
saveas(gcf,[figureDir 'percentCompletedPR_zoom_driftRL_valueUpdate.fig'],'fig')
print([figureDir 'percentCompletedPR_zoom_driftRL_valueUpdate.png'],'-dpng','-r600')

%% Example drift diffusion
oneHitBound = 0;
noiseAmplitude = .1;
req = 60;
driftRateCoeff = .01;
Qalt = .25;
Qcur = .4;
count = 0;
while (~oneHitBound && count <= 100)
    [trajectory1,hitBound1,boundcross1,Qs1] = drift_valueUpdateProcess(Qcur,Qalt,driftRateCoeff,req,noiseAmplitude,0,1,60);
    [trajectory2,hitBound2,boundcross2,Qs2] = drift_valueUpdateProcess(Qcur,Qalt,driftRateCoeff,req,noiseAmplitude,0,1,60);
    if (hitBound1 || hitBound2 & ~(hitBound1 & hitBound2))
        figure;
        plot(1:req,trajectory1,'b','linewidth',3)
        hold on;
        plot(1:req,trajectory2,'r','linewidth',3)
        ylabel('X(t)','fontsize',15,'fontweight','bold')
        xlabel('Lever press','fontsize',15,'fontweight','bold')
        ylim([min([min(trajectory1) min(trajectory2)]),1.2])
        plot(1:req,ones(1,req),'--')
        set(gcf,'Position',[10 10 1600 1200])
        saveas(gcf,[figureDir 'driftDiffusion_example.fig'],'fig')
        print([figureDir 'driftDiffusion_example.png'],'-dpng','-r600')
        oneHitBound = true;
    end
    count=count+1;
end

%% BIC for driftRL_valueUpdate and logistic models
BICs = [];
nDataPoints = 160;
for i=1:4
    [scores,allParams] = getScores(driftRL_valueUpdateDirs{i},'useSave',true);
    [~,minInd] = min(scores);
    params = load([driftRL_valueUpdateDirs{1} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end

for i=1:4
    [scores,allParams] = getScores(logisticRLDirs{i},'useSave',true);
    [~,minInd] = min(scores);
    params = load([logisticRLDirs{1} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end
xtick = [1:4 6:9];
xticklabels = {'DDM',{'DDM','no fatigue'},{'DDM','no WF'},{'DDM','no fatigue','no WF'},...
    'logistic',{'logistic','no fatigue'},{'logistic','no WF'},{'logistic','no fatigue','no WF'}};
figure;
bar(xtick,BICs)
set(gca,'xtick',xtick,'fontsize',12,'fontweight','bold')
set(gcf,'Position',[10 10 1400 800])
ylabel('BIC','fontsize',20,'fontweight','bold')
saveas(gcf,[figureDir '/model_comparison_BIC.fig'],'fig')
print([figureDir '/model_comparison_BIC_driftRL_valueUpdate_logistic.png'],'-dpng','-r600')

%% BIC for driftRL and logistic models
BICs = [];
nDataPoints = 160;
for i=1:4
    [scores,allParams] = getScores(driftRLDirs{i},'useSave',true);
    [~,minInd] = min(scores);
    params = load([driftRLDirs{1} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end

for i=1:4
    [scores,allParams] = getScores(logisticRLDirs{i},'useSave',true);
    [~,minInd] = min(scores);
    params = load([logisticRLDirs{1} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end
xtick = [1:4 6:9];
xticklabels = {'DDM',{'DDM','no fatigue'},{'DDM','no WF'},{'DDM','no fatigue','no WF'},...
    'logistic',{'logistic','no fatigue'},{'logistic','no WF'},{'logistic','no fatigue','no WF'}};
figure;
bar(xtick,BICs)
set(gca,'xtick',xtick,'fontsize',12,'fontweight','bold')
set(gcf,'Position',[10 10 1400 800])
ylabel('BIC','fontsize',20,'fontweight','bold')
saveas(gcf,[figureDir '/model_comparison_BIC.fig'],'fig')
print([figureDir '/model_comparison_BIC_driftRL_logistic.png'],'-dpng','-r600')