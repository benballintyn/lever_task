%pizza_talk_figs_2021
%%
figureDir = '/home/ben/phd/lever_task/publication/modeling/12_7_21';
if (~exist(figureDir,'dir'))
    mkdir(figureDir)
end

useOnly120Trials = true;

% Check which computer is in use --> which external drive to use
[status,result] = system('hostname');
computerName = strtrim(result);
if (strcmp(computerName,'miller-lab-ubuntu2'))
    externalHDDir = '/media/ben/Manwe/';
elseif (strcmp(computerName,'silmaril'))
    externalHDDir = '/media/ben/Varda/';
else
    error('hostname not recognized')
end

externalDataDir = [externalHDDir 'phd/lever_task/'];
if (useOnly120Trials)
    driftRLDir = [externalDataDir 'driftRL/results/Only120Trials/logprob_joint/value_based_drift/'];
    driftSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues';
else
    driftRLDir = [externalDataDir 'driftRL/results/logprob_joint/value_based_drift/'];
    driftSimDataDirs{1} = 'bandit_e_greedy_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{2} = 'bandit_e_greedy_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{3} = 'bandit_e_greedy_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{4} = 'bandit_e_greedy_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{5} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
end

for i=1:length(driftSimDataDirs)
    driftRLDirs{i} = [driftRLDir driftSimDataDirs{i} '/optimized/'];
end

if (useOnly120Trials)
    logisticRLDir = [externalDataDir 'logisticAbortRL/results/Only120Trials/logprob_joint/'];
    logisticSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{5} = 'bandit_e_greedy_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{6} = 'bandit_e_greedy_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{7} = 'bandit_e_greedy_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{8} = 'bandit_e_greedy_initialization_mean_reward_forgettingType_decayToInitialValues';
else
    logisticRLDir = [externalDataDir 'logisticAbortRL/results/logprob_joint/'];
    logisticSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{5} = 'bandit_e_greedy_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
end

for i=1:length(logisticSimDataDirs)
    logisticRLDirs{i} = [logisticRLDir logisticSimDataDirs{i} '/optimized/'];
end

if (useOnly120Trials)
    driftRL_valueUpdateDir = [externalDataDir 'driftRL_valueUpdate/results/Only120Trials/logprob_joint/value_based_drift/'];
    driftRL_valueUpdateSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftRL_valueUpdateSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftRL_valueUpdateSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftRL_valueUpdateSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues';
else
    driftRL_valueUpdateDir = [externalDataDir 'driftRL_valueUpdate/results/logprob_joint/value_based_drift/'];
    driftRL_valueUpdateSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftRL_valueUpdateSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftRL_valueUpdateSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftRL_valueUpdateSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues';
end

for i=1:length(driftRL_valueUpdateSimDataDirs)
    driftRL_valueUpdateDirs{i} = [driftRL_valueUpdateDir driftRL_valueUpdateSimDataDirs{i} '/optimized/'];
end

%% DriftRL figs first
for i=1:length(driftRLDirs)
    disp(['Making plots for: ' driftRLDirs{i}])
    [scores,allParams] = plotResults(driftRLDirs{i},1,'fitStatsOnly',true,'Only120Trials',useOnly120Trials);
    saveas(gcf,[figureDir 'driftRL_' driftSimDataDirs{i} '_bestFitStats.fig'],'fig')
    print([figureDir 'driftRL_' driftSimDataDirs{i} '_bestFitStats.png'],'-dpng','-r600')
end

%% LogisticAbort figs
for i=1:length(logisticRLDirs)
    disp(['Making plots for: ' logisticRLDirs{i}])
    [scores,allParams] = plotResults(logisticRLDirs{i},1,'fitStatsOnly',true,'Only120Trials',useOnly120Trials);
    saveas(gcf,[figureDir 'logisticAbortRL_' logisticSimDataDirs{i} '_bestFitStats.fig'],'fig')
    print([figureDir 'logisticAbortRL_' logisticSimDataDirs{i} '_bestFitStats.png'],'-dpng','-r600')
end

%% driftRL_valueUpdate figs
for i=1:length(driftRL_valueUpdateDirs)
    disp(['Making plots for: ' driftRL_valueUpdateDirs{i}])
    [scores,allParams] = plotResults(driftRL_valueUpdateDirs{i},1,'fitStatsOnly',true,'Only120Trials',useOnly120Trials);
    saveas(gcf,[figureDir 'driftRL_valueUpdate_' driftRL_valueUpdateSimDataDirs{i} '_bestFitStats.fig'],'fig');
    print([figureDir 'driftRL_valueUpdate_' driftRL_valueUpdateSimDataDirs{i} '_bestFitStats.png'],'-dpng','-r600')
end

%% Aborted trials analysis
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
nTrialsAborted = load('optimality/WT/nTrialsAborted.mat'); nTrialsAborted=nTrialsAborted.nTrialsAborted;
nTrialsCompleted = load('optimality/WT/nTrialsCompleted.mat'); nTrialsCompleted=nTrialsCompleted.nTrialsCompleted;
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
% Logistic RL first (model without fatigue)
%[scores,allParams] = plotResults(logisticRLDirs{1},1,'fitStatsOnly',true);
[scores,allParams] = getScores(logisticRLDirs{2});
[~,maxInd] = min(scores);
percentCompletedPR = load([logisticRLDirs{2} '/' num2str(maxInd) '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
figure;
for i=1:4
    subplot(2,2,i)
    curvals = percentCompletedPR_allTrials{i}(percentCompletedPR_allTrials{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    hold on;
    curvals = percentCompletedPR{i}(percentCompletedPR{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    title(sessionTypes{i},'fontsize',15,'fontweight','bold')
    xlabel('Fraction of PR trial completed','fontsize',15,'fontweight','bold')
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

% Logistic model no WF
[scores,allParams] = getScores(logisticRLDirs{4});
[~,maxInd] = min(scores);
percentCompletedPR = load([logisticRLDirs{4} '/' num2str(maxInd) '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
figure;
for i=1:4
    subplot(2,2,i)
    curvals = percentCompletedPR_allTrials{i}(percentCompletedPR_allTrials{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    hold on;
    curvals = percentCompletedPR{i}(percentCompletedPR{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    title(sessionTypes{i},'fontsize',15,'fontweight','bold')
    xlabel('Fraction of PR trial completed','fontsize',15,'fontweight','bold')
    if (i == 1)
        legend({'Mice','Model'},'fontsize',15,'fontweight','bold','Location','northeast')
    end
    if (i == 1 || i == 3)
        ylabel('Probability density','fontsize',15,'fontweight','bold')
    end
end
set(gcf,'Position',[10 10 1600 1600])
saveas(gcf,[figureDir 'percentCompletedPR_zoom_logisticRL_noWF.fig'],'fig')
print([figureDir 'percentCompletedPR_zoom_logisticRL_noWF.png'],'-dpng','-r600')

% Drift RL next (model without fatigue)
%[scores,allParams] = plotResults(driftRLDirs{1},1,'fitStatsOnly',true);
[scores,allParams] = getScores(driftRLDirs{2});
[~,maxInd] = min(scores);
percentCompletedPR = load([driftRLDirs{2} '/' num2str(maxInd) '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
figure;
for i=1:4
    subplot(2,2,i)
    curvals = percentCompletedPR_allTrials{i}(percentCompletedPR_allTrials{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    hold on;
    curvals = percentCompletedPR{i}(percentCompletedPR{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    title(sessionTypes{i},'fontsize',15,'fontweight','bold')
    xlabel('Fraction of PR trial completed','fontsize',15,'fontweight','bold')
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

% Drift RL model w/o WF
[scores,allParams] = getScores(driftRLDirs{4});
[~,maxInd] = min(scores);
percentCompletedPR = load([driftRLDirs{4} '/' num2str(maxInd) '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
figure;
for i=1:4
    subplot(2,2,i)
    curvals = percentCompletedPR_allTrials{i}(percentCompletedPR_allTrials{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    hold on;
    curvals = percentCompletedPR{i}(percentCompletedPR{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    title(sessionTypes{i},'fontsize',15,'fontweight','bold')
    xlabel('Fraction of PR trial completed','fontsize',15,'fontweight','bold')
    if (i == 1)
        legend({'Mice','Model'},'fontsize',15,'fontweight','bold','Location','northeast')
    end
    if (i == 1 || i == 3)
        ylabel('Probability density','fontsize',15,'fontweight','bold')
    end
end
set(gcf,'Position',[10 10 1600 1600])
saveas(gcf,[figureDir 'percentCompletedPR_zoom_driftRL_noWF.fig'],'fig')
print([figureDir 'percentCompletedPR_zoom_driftRL_noWF.png'],'-dpng','-r600')

% DriftRL_valueUpdate last (model without fatigue)
%[scores,allParams] = plotResults(driftRL_valueUpdateDirs{1},1,'fitStatsOnly',true);
[scores,allParams] = getScores(driftRL_valueUpdateDirs{2});
[~,maxInd] = min(scores);
percentCompletedPR = load([driftRL_valueUpdateDirs{2} '/' num2str(maxInd) '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
figure;
for i=1:4
    subplot(2,2,i)
    curvals = percentCompletedPR_allTrials{i}(percentCompletedPR_allTrials{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    hold on;
    curvals = percentCompletedPR{i}(percentCompletedPR{i} < 1);
    histogram(curvals,'normalization','pdf','binwidth',.02)
    title(sessionTypes{i},'fontsize',15,'fontweight','bold')
    xlabel('Fraction of PR trial completed','fontsize',15,'fontweight','bold')
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
    params = load([driftRL_valueUpdateDirs{i} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end

for i=1:4
    [scores,allParams] = getScores(logisticRLDirs{i},'useSave',true);
    [~,minInd] = min(scores);
    params = load([logisticRLDirs{i} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end
xtick = [1:4 6:9];
xticklabels = {'DDM_{VU}','DDM_{VU} F-','DDM_{VU} WF-','DDM_{VU} F-/WF-',...
    'LOG(\tau)','LOG(\tau) F-','LOG(\tau) WF-','LOG(\tau) F-/WF-'};
figure;
bar(xtick,BICs)
set(gca,'xtick',xtick,'xticklabels',xticklabels,'fontsize',12,'fontweight','bold')
set(gcf,'Position',[10 10 1800 800])
ylabel('BIC','fontsize',20,'fontweight','bold')
saveas(gcf,[figureDir '/model_comparison_BIC_driftRL_valueUpdate_logistic.fig'],'fig')
print([figureDir '/model_comparison_BIC_driftRL_valueUpdate_logistic.png'],'-dpng','-r600')

%% BIC for driftRL and logistic models
BICs = [];
nDataPoints = 160;
for i=1:4
    [scores,allParams] = getScores(driftRLDirs{i},'useSave',true);
    [~,minInd] = min(scores);
    params = load([driftRLDirs{i} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end

for i=1:4
    [scores,allParams] = getScores(logisticRLDirs{i},'useSave',true);
    [~,minInd] = min(scores);
    params = load([logisticRLDirs{i} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end
xtick = [1:4 6:9];
xticklabels = {'DDM','DDM F-','DDM WF-','DDM F-/WF-',...
    'LOG(\tau)','LOG(\tau) F-','LOG(\tau) WF-','LOG(\tau) F-/WF-'};
figure;
bar(xtick,BICs)
set(gca,'xtick',xtick,'xticklabels',xticklabels,'fontsize',12,'fontweight','bold')
set(gcf,'Position',[10 10 1800 800])
ylabel('BIC','fontsize',20,'fontweight','bold')
saveas(gcf,[figureDir '/model_comparison_BIC_driftRL_logistic.fig'],'fig')
print([figureDir '/model_comparison_BIC_driftRL_logistic.png'],'-dpng','-r600')

%% BIC for logisticRL softmax and e-greedy
BICs = [];
nDataPoints = 160;
for i=1:8
    [scores,allParams] = getScores(logisticRLDirs{i},'useSave',true);
    [~,minInd] = min(scores);
    params = load([logisticRLDirs{i} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end

xtick = [1:4 6:9];
xticklabels = {'LOG(\tau)','LOG(\tau) F-','LOG(\tau) WF-','LOG(\tau) F-/WF-',...
    'LOG(\epsilon)','LOG(\epsilon) F-','LOG(\epsilon) WF-','LOG(\epsilon) F-/WF-'};
figure;
bar(xtick,BICs)
set(gca,'xtick',xtick,'xticklabels',xticklabels,'fontsize',12,'fontweight','bold')
set(gcf,'Position',[10 10 1800 800])
ylabel('BIC','fontsize',20,'fontweight','bold')
saveas(gcf,[figureDir '/model_comparison_BIC_logistic.fig'],'fig')
print([figureDir '/model_comparison_BIC_logistic.png'],'-dpng','-r600')

%% BIC for all drift models and F- models for driftRL_valueUpdate and logistic models
BICs = [];
nDataPoints = 160;
for i=1:4
    [scores,allParams] = getScores(driftRLDirs{i},'useSave',true);
    [~,minInd] = min(scores);
    params = load([driftRLDirs{i} '/' num2str(minInd) '/params.mat']); params=params.params;
    curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
    BICs = [BICs curBIC];
end

[scores,allParams] = getScores(logisticRLDirs{2},'useSave',true);
[~,minInd] = min(scores);
params = load([logisticRLDirs{i} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

[scores,allParams] = getScores(logisticRLDirs{6},'useSave',true);
[~,minInd] = min(scores);
params = load([logisticRLDirs{i} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

[scores,allParams] = getScores(driftRL_valueUpdateDirs{2},'useSave',true);
[~,minInd] = min(scores);
params = load([driftRL_valueUpdateDirs{i} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

xtick = [1:4 6:7 9];
xticklabels = {'DDM','DDM F-','DDM WF-','DDM F-/WF-',...
    'LOG(\tau) F-','LOG(\epsilon) F-','DDM_{VU} F-'};
figure;
bar(xtick,BICs)
set(gca,'xtick',xtick,'xticklabels',xticklabels,'fontsize',12,'fontweight','bold')
set(gcf,'Position',[10 10 1800 800])
ylabel('BIC','fontsize',20,'fontweight','bold')
saveas(gcf,[figureDir '/model_comparison_BIC_decent_models.fig'],'fig')
print([figureDir '/model_comparison_BIC_decent_models.png'],'-dpng','-r600')