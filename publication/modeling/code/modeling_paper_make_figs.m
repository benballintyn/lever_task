%% MODELING PAPER FIGURES
SAVE_DIR = '~/phd/lever_task/publication/modeling/figures/v1/';
if (~exist(SAVE_DIR,'dir'))
    mkdir(SAVE_DIR)
end

NL_optimal = [10 22 28 58];
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
NL_observed = load(['optimality/WT/NL_observed.mat']); NL_observed=NL_observed.NL_observed;

useOnly120Trials = true;

externalHDDir = '/media/ben/Varda/';
externalDataDir = [externalHDDir 'phd/lever_task/'];
if (useOnly120Trials)
    driftRLDir = [externalDataDir 'driftRL/results/Only120Trials/logprob_joint/value_based_drift/'];
    driftSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    driftSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues';
else
    driftRLDir = [externalDataDir 'driftRL/results/logprob_joint/value_based_drift/'];
    driftSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
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
else
    logisticRLDir = [externalDataDir 'logisticAbortRL/results/logprob_joint/'];
    logisticSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
    logisticSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues';
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


%% Fig 1C (Optimal switching points with FR/PR values vs. trial #)
load('~/phd/lever_task/publication/modeling/data_files/perfectAgentQ.mat');
for i=1:4
    subplot(2,2,i)
    plot(Q{i})
    hold on;
    plot([NL_optimal(i) NL_optimal(i)],[0 8],'Color','k')
    legend({'FR','PR','Optimal PR #'})
    title(sessionTypes{i})
    ylim([0 8])
    xlim([0 60])
    xlabel('PR trial #','fontsize',20,'fontweight','bold')
    ylabel('EoR','fontsize',20,'fontweight','bold')
    set(gca,'fontsize',15,'fontweight','bold')
end
set(gcf,'Position',[10 10 1600 1200])
figName = 'perfect_agent_Qvalues_optimal_labeled';
saveas(gcf,[SAVE_DIR figName '.fig'],'fig')
print([SAVE_DIR figName '.png'],'-dpng','-r600')

%% Fig 2A (P(FR) & P(PR) for mouse and best fit agent)

% Mouse data first
pS = load('~/phd/lever_task/optimality/WT/pS.mat'); pS=pS.pS;
pL = load('~/phd/lever_task/optimality/WT/pL.mat'); pL=pL.pL;
figure;
for i=1:4
    subplot(1,4,i)
    plot(pS{i}(1:120)); hold on; plot(pL{i}(1:120))
    xlabel('Trial number','fontsize',20,'fontweight','bold')
    ylabel('Choice probability','fontsize',20,'fontweight','bold')
    legend({'P(FR)','P(PR)'},'location','best','fontsize',15,'fontweight','bold')
    title(sessionTypes{i})
    
    % Expand axes to eliminate whitespace
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end
set(gcf,'Position',[10 10 2000 500])
figName = 'p_FR_p_PR_mouse';
saveas(gcf,[SAVE_DIR figName '.fig'],'fig')
print([SAVE_DIR figName '.png'],'-dpng','-r600')

% Best logisticRL agent
[scores,allParams] = getScores(logisticRLDirs{2},'useSave',true);
[~,minInd] = min(scores);
pS = load([logisticRLDirs{2} '/' num2str(minInd) '/pS.mat']); pS=pS.pS;
pL = load([logisticRLDirs{2} '/' num2str(minInd) '/pL.mat']); pL=pL.pL;
figure;
for i=1:4
    subplot(1,4,i)
    plot(pS{i}(1:120)); hold on; plot(pL{i}(1:120))
    xlabel('Trial number','fontsize',20,'fontweight','bold')
    ylabel('Choice probability','fontsize',20,'fontweight','bold')
    legend({'P(FR)','P(PR)'},'location','best','fontsize',15,'fontweight','bold')
    title(sessionTypes{i})
    
    % Expand axes to eliminate whitespace
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end
set(gcf,'Position',[10 10 2000 500])
figName = 'p_FR_p_PR_agent';
saveas(gcf,[SAVE_DIR figName '.fig'],'fig')
print([SAVE_DIR figName '.png'],'-dpng','-r600')

%%