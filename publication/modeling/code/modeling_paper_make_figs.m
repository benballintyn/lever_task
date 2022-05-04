%% MODELING PAPER FIGURES
clear all;
BASE_DIR = '~/phd/lever_task/publication/modeling/';
FIGURE_DIR = [BASE_DIR 'figures/thesis/'];
DATA_DIR = [BASE_DIR 'data_files/'];
if (~exist(BASE_DIR,'dir'))
    mkdir(BASE_DIR)
end
if (~exist(FIGURE_DIR,'dir'))
    mkdir(FIGURE_DIR)
end
if (~exist(DATA_DIR,'dir'))
    mkdir(DATA_DIR)
end

NL_optimal = [10 22 28 58];
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
NL_observed = load(['optimality/WT/NL_observed.mat']); NL_observed=NL_observed.NL_observed;

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
    driftSimDataDirs{5} = 'bandit_UCB_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues_noAbortANS';
else
    driftRLDir = [externalDataDir 'driftRL/results/logprob_joint/value_based_drift/'];
    driftSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues';
end

for i=1:length(driftSimDataDirs)
    driftRLDirs{i} = [driftRLDir driftSimDataDirs{i} '/optimized/'];
end

if (useOnly120Trials)
    logisticRLDir = [externalDataDir 'logisticAbortRL/results/Only120Trials/logprob_joint/'];
    logisticSimDataDirs{1} = 'bandit_softmax_ansUtilityFunc_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues_noAbortANS';
    logisticSimDataDirs{2} = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues_noAbortANS';
    logisticSimDataDirs{3} = 'bandit_softmax_pressUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues_noAbortANS';
    logisticSimDataDirs{4} = 'bandit_softmax_initialization_mean_reward_forgettingType_decayToInitialValues_noAbortANS';
    logisticSimDataDirs{5} = 'bandit_UCB_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues_noAbortANS';
    logisticSimDataDirs{6} = 'bandit_e_greedy_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues_noAbortANS';
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
figure;
elim_whitespace;
for i=1:4
    subplot(2,2,i)
    plot(Q{i})
    hold on;
    plot([NL_optimal(i) NL_optimal(i)],[0 8],'Color','k')
    if (i == 1)
        legend({'FR','PR','Optimal PR #'})
    end
    title(sessionTypes{i})
    ylim([0 8])
    xlim([0 60])
    xlabel('PR trial #','fontsize',20,'fontweight','bold')
    ylabel('EoR','fontsize',20,'fontweight','bold')
    set(gca,'fontsize',15,'fontweight','bold')
end
set(gcf,'Position',[10 10 1600 1200])
figName = 'fig1C_perfect_agent_Qvalues_optimal_labeled';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig 2A (P(FR) & P(PR) for mouse and best fit agent)

% Mouse data first
pS = load('~/phd/lever_task/optimality/WT/pS.mat'); pS=pS.pS;
pL = load('~/phd/lever_task/optimality/WT/pL.mat'); pL=pL.pL;
figure;
elim_whitespace
for i=1:4
    subplot(1,4,i)
    plot(pS{i}(1:120)); 
    hold on; 
    plot(pL{i}(1:120))
    plot([NL_optimal(i) NL_optimal(i)],[0 1],'k','linewidth',2)
    xlabel('Trial number','fontsize',20,'fontweight','bold')
    ylabel('Choice probability','fontsize',20,'fontweight','bold')
    if (i == 1)
        legend({'P(FR)','P(PR)'},'location','east','fontsize',15,'fontweight','bold')
    end
    title(sessionTypes{i})
end
set(gcf,'Position',[10 10 2000 500])
figName = 'fig2A_p_FR_p_PR_mouse';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

% Best logisticRL agent
[scores,allParams] = getScores(logisticRLDirs{2});
[~,minInd] = min(scores);
pS = load([logisticRLDirs{2} '/' num2str(minInd) '/pS.mat']); pS=pS.pS;
pL = load([logisticRLDirs{2} '/' num2str(minInd) '/pL.mat']); pL=pL.pL;
figure;
elim_whitespace;
for i=1:4
    subplot(1,4,i)
    plot(pS{i}(1:120));
    hold on;
    plot(pL{i}(1:120))
    plot([NL_optimal(i) NL_optimal(i)],[0 1],'k','linewidth',2)
    
    xlabel('Trial number','fontsize',20,'fontweight','bold')
    ylabel('Choice probability','fontsize',20,'fontweight','bold')
    if (i == 1)
        legend({'P(FR)','P(PR)'},'location','east','fontsize',15,'fontweight','bold')
    end
    title(sessionTypes{i})
end
set(gcf,'Position',[10 10 2000 500])
figName = 'fig2B_p_FR_p_PR_agent';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig 3A (Fraction trials aborted vs. relative value w/ sigmoidal fit)
% Load mouse data and do abort analysis
mice = load('processed_data/WT/animals.mat'); mice=mice.animals;
abort_countPR_all = zeros(100,100);
total_countPR_all = zeros(100,100);
abort_countFR_all = zeros(2,12);
total_countFR_all = zeros(2,12);

abort_countPR = cell(1,4);
total_countPR = cell(1,4);
abort_countFR = cell(1,4);
total_countFR = cell(1,4);

percentCompletedPR = cell(1,4);
percentCompletedFR = cell(1,4);
nAborted = cell(1,4);
dataMatrix = [];
for i=1:4
    abort_countPR{i} = zeros(100,100);
    total_countPR{i} = zeros(100,100);
    if (i == 1 || i == 3)
        abort_countFR{i} = zeros(1,6);
        total_countFR{i} = zeros(1,6);
    elseif (i == 2 || i == 4)
        abort_countFR{i} = zeros(1,12);
        total_countFR{i} = zeros(1,12);
    end
end

count = 0;
mouseNames = {};
for i=1:length(mice)
    data = load(['processed_data/WT/' mice{i} '_ReProcessedData.mat']); data=data.ProcessedData;
    for j=1:length(data)
        sessType = [num2str(data{j}.Reward_L / data{j}.Reward_S) 'xFR' num2str(data{j}.NumPressRequired_S)];
        switch sessType
            case '2xFR6'
                sessInd = 1;
                FRind = 1;
            case '2xFR12'
                sessInd = 2;
                FRind = 2;
            case '5xFR6'
                sessInd = 3;
                FRind = 1;
            case '5xFR12'
                sessInd = 4;
                FRind = 2;
        end
        dayAbortedCount = 0;
        for k=1:data{j}.TotalTrialsCompleted
            count = count+1;
            mouseNames{count} = data{j}.AnimalName;
            relativeValue = (data{j}.Reward_L/data{j}.NumPressRequired_L(k)) - (data{j}.Reward_S/data{j}.NumPressRequired_S);
            if (isnan(data{j}.SideRewarded(k)))
                dayAbortedCount = dayAbortedCount + 1;
                if (data{j}.SideChosen(k) == 0)
                    abortPress = data{j}.LeverPressesByTrial(k); % - 1 to take out press that represents side choice
                    pressRequirement = data{j}.NumPressRequired_S; % - 1 to take out press that represents side choice
                    percentCompletedFR{sessInd} = [percentCompletedFR{sessInd} abortPress/pressRequirement];
                    
                    abort_countFR_all(FRind,abortPress) = abort_countFR_all(FRind,abortPress) + 1;
                    total_countFR_all(FRind,1:abortPress) = total_countFR_all(FRind,1:abortPress) + 1;
                    
                    abort_countFR{sessInd}(abortPress) = abort_countFR{sessInd}(abortPress) + 1;
                    total_countFR{sessInd}(1:abortPress) = total_countFR{sessInd}(1:abortPress) + 1;
                    
                    dataMatrix(:,count) = [data{j}.SideChosen(k); isnan(data{j}.SideRewarded(k)); data{j}.LeverPressesByTrial(k); abortPress/pressRequirement; relativeValue; sessInd; pressRequirement; k; j; data{j}.Reward_L/data{j}.NumPressRequired_L(k); data{j}.Reward_S/data{j}.NumPressRequired_S];
                elseif (data{j}.SideChosen(k) == 1)
                    abortPress = data{j}.LeverPressesByTrial(k);
                    pressRequirement = data{j}.NumPressRequired_L(k);
                    percentCompletedPR{sessInd} = [percentCompletedPR{sessInd} abortPress/pressRequirement];
                    
                    abort_countPR_all(pressRequirement,abortPress) = abort_countPR_all(pressRequirement,abortPress) + 1;
                    total_countPR_all(pressRequirement,1:abortPress) = total_countPR_all(pressRequirement,1:abortPress) + 1;
                    
                    abort_countPR{sessInd}(pressRequirement,abortPress) = abort_countPR{sessInd}(pressRequirement,abortPress) + 1;
                    total_countPR{sessInd}(pressRequirement,1:abortPress) = total_countPR{sessInd}(pressRequirement,1:abortPress) + 1;
                    
                    dataMatrix(:,count) = [data{j}.SideChosen(k); isnan(data{j}.SideRewarded(k)); data{j}.LeverPressesByTrial(k); abortPress/pressRequirement; relativeValue; sessInd; pressRequirement; k; j; data{j}.Reward_L/data{j}.NumPressRequired_L(k); data{j}.Reward_S/data{j}.NumPressRequired_S];
                elseif (isnan(data{j}.SideChosen(k)))
                    continue;
                else
                    error(['mouse ' mice{i} ' data{' num2str(j) '}.SideChosen(' num2str(k) ') is not valid'])
                end
            else
                if (data{j}.SideChosen(k) == 0)
                    pressRequirement = data{j}.NumPressRequired_S;
                    
                    total_countFR_all(FRind,1:pressRequirement) = total_countFR_all(FRind,1:pressRequirement) + 1;
                    total_countFR{sessInd} = total_countFR{sessInd} + 1;
                    
                    dataMatrix(:,count) = [data{j}.SideChosen(k); isnan(data{j}.SideRewarded(k)); data{j}.LeverPressesByTrial(k); 1; relativeValue; sessInd; pressRequirement; k; j; data{j}.Reward_L/data{j}.NumPressRequired_L(k); data{j}.Reward_S/data{j}.NumPressRequired_S];
                elseif (data{j}.SideChosen(k) == 1)
                    pressRequirement = data{j}.NumPressRequired_L(k);
                    
                    total_countPR_all(pressRequirement,1:pressRequirement) = total_countPR_all(pressRequirement,1:pressRequirement) + 1;
                    total_countPR{sessInd}(pressRequirement,1:pressRequirement) = total_countPR{sessInd}(pressRequirement,1:pressRequirement) + 1;
                    
                    dataMatrix(:,count) = [data{j}.SideChosen(k); isnan(data{j}.SideRewarded(k)); data{j}.LeverPressesByTrial(k); 1; relativeValue; sessInd; pressRequirement; k; j; data{j}.Reward_L/data{j}.NumPressRequired_L(k); data{j}.Reward_S/data{j}.NumPressRequired_S];
                elseif (isnan(data{j}.SideChosen(k)))
                    continue;
                else
                    error(['mouse ' mice{i} ' data{' num2str(j) '}.SideChosen(' num2str(k) ') is not valid'])
                end
            end
        end
        nAborted{sessInd} = [nAborted{sessInd} dayAbortedCount];
    end
end

dataTable = array2table(dataMatrix','VariableNames',{'SideChosen','WasAborted','LeverPressesByTrial','FracCompleted','RelativeValue','SessionIndex','PressRequirment','TrialNum','DayIndex', 'PRvalue', 'FRvalue'});
dataTable.mouseNames = mouseNames';

Pabort_PR_all = abort_countPR_all./total_countPR_all;
Pabort_FR_all = abort_countFR_all./total_countFR_all;
for i=1:4
    Pabort_PR{i} = abort_countPR{i}./total_countPR{i};
    Pabort_FR{i} = abort_countFR{i}./total_countFR{i};
end

save([DATA_DIR 'percentCompletedFR.mat'],'percentCompletedFR','-mat')
save([DATA_DIR 'percentCompletedPR.mat'],'percentCompletedPR','-mat')
save([DATA_DIR 'Pabort_FR.mat'],'Pabort_FR','-mat')
save([DATA_DIR 'Pabort_PR.mat'],'Pabort_PR','-mat')
save([DATA_DIR 'Pabort_PR_all.mat'],'Pabort_PR_all','-mat')
save([DATA_DIR 'Pabort_FR_all.mat'],'Pabort_FR_all','-mat')
save([DATA_DIR 'dataMatrix.mat'],'dataMatrix','-mat')

bins = [-.4:.1:4.9; -.3:.1:5];
nbins = length(bins);
meanbins = mean(bins,1);
PRinds = find(dataMatrix(1,:) == 1);
fracAborted = zeros(1,nbins);
ninds = zeros(1,nbins);
for i=1:size(bins,2)
    inds = find(dataMatrix(5,:) > bins(1,i) & dataMatrix(5,:) <= bins(2,i));
    prinds = intersect(PRinds,inds);
    fracAborted(i) = sum(dataMatrix(2,prinds))/length(prinds);
end

PRtrials = dataTable.SideChosen == 1;
NaNinds = isnan(fracAborted);
%ft = fittype(' - 1/(1 + exp(-A*x)) + B');
%F = fit(meanbins(~NaNinds)',fracAborted(~NaNinds)',ft);
f = @(beta,x) 1 - beta(1) - (beta(2) - beta(1))./(beta(3) + exp(-beta(4)*(x - beta(5))));
beta0 = [0 1 1 1 0];
%NLM = fitnlm(meanbins(~NaNinds),fracAborted(~NaNinds),f,beta0);
NLM = fitnlm(dataTable.RelativeValue(PRtrials),dataTable.WasAborted(PRtrials),f,beta0);
%Fall = fit(dataTable.RelativeValue(PRtrials),dataTable.WasAborted(PRtrials),ft);
%f = @(x,A,B) 1 - 1./(1 + exp(-A*x)) + B;
x = -1:.001:5;
figure;
%elim_whitespace
scatter(meanbins,fracAborted,100,'k.'); 
hold on;
%plot(x,f(x,F.A,F.B))
%plot(x,f(x,Fall.A,Fall.B))
plot(x,NLM.predict(x'))

legend({'Binned data',...
       ['logistic']},'fontsize',20,'fontweight','bold')
xlabel('Relative value (PR - FR)','FontSize',20,'FontWeight','bold')
ylabel('Fraction of trials aborted','FontSize',20,'FontWeight','bold')
set(gcf,'Position',[10 10 1600 1200])
figName = 'fig3A_fracAborted_vs_relativeValue';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig 3B (P(abort) after % of trial completed (depends on prior section)
for i=1:4
    sessInds = find(dataTable.SessionIndex == i);
    prInds = find(dataTable.SideChosen == 1);
    curInds = intersect(sessInds,prInds);
    [n_allTrials,~] = histcounts(dataTable.FracCompleted(curInds),100);
    n_allTrials(end) = 0;
    percentCompletedCDFs_PR_allTrials{i} = cumsum(n_allTrials)/length(curInds);
    pdfs_PR_allTrials{i} = n_allTrials./length(curInds);
    survivor_PR_allTrials{i} = 1 - percentCompletedCDFs_PR_allTrials{i};
end
for i=1:4
    for j=1:(length(survivor_PR_allTrials{i}) - 1)
        abortConditionalProb_PR_allTrials{i}(j) = (survivor_PR_allTrials{i}(j) - survivor_PR_allTrials{i}(j+1))/survivor_PR_allTrials{i}(j);
    end
end

figure;
elim_whitespace;
for i=1:4
    subplot(2,2,i)
    plot(abortConditionalProb_PR_allTrials{i})
    xlabel('% of trial completed','fontsize',15,'fontweight','bold')
    ylabel({'Conditional probability of','abort after x% of trial'},'fontsize',15,'fontweight','bold')
end
set(gcf,'Position',[10 10 1400 1200])
figName = 'fig3B_PR_trial_conditional_probability_of_abort_after_press_x_allTrials';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig 3C Drift Diffusion example
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
    if ((hitBound1 || hitBound2) & ~(hitBound1 & hitBound2))
        figure;
        elim_whitespace
        plot(1:req,trajectory1,'b','linewidth',3)
        hold on;
        plot(1:req,trajectory2,'r','linewidth',3)
        ylabel('Decision variable','fontsize',20,'fontweight','bold')
        xlabel('Lever press #','fontsize',20,'fontweight','bold')
        ylim([min([min(trajectory1) min(trajectory2)]),1.2])
        plot(1:req,ones(1,req),'--')
        set(gcf,'Position',[10 10 1600 1200])
        figName = 'fig3C_driftDiffusion_example';
        saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
        saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
        print([FIGURE_DIR figName '.png'],'-dpng','-r600')
        oneHitBound = true;
    end
    count=count+1;
end

%% Fig 4A % PR trial completed for mice and model - DDM
[scores,allParams] = plotResults(driftRLDirs{2},1,'Only120Trials',true,'PercentAbortedOnly',true);
figName = 'fig4A_driftRL_percent_PR_completed';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig 4B % PR trial completed for mice and model - Logistic
[scores,allParams] = plotResults(logisticRLDirs{2},1,'Only120Trials',true,'PercentAbortedOnly',true);
figName = 'fig4B_logisticRL_percent_PR_completed';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig 5A Best fit drift RL stats with Weber-Fechner
[scores,allParams] = plotResults(driftRLDirs{2},1,'fitStatsOnly',true,'Only120Trials',true);
figName = 'fig5A_driftRL_WF_bestFit';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig5B Best fit driftRL stats without Weber-Fechner
[scores,allParams] = plotResults(driftRLDirs{4},1,'fitStatsOnly',true,'Only120Trials',true);
figName = 'fig5B_driftRL_noWF_bestFit';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig 6A Best fit logisticAbortRL stats with Weber-Fechner
[scores,allParams] = plotResults(logisticRLDirs{2},1,'fitStatsOnly',true,'Only120Trials',true);
figName = 'fig6A_logisticAbortRL_WF_bestFit';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig 6B Best fit logisticAbortRL stats without Weber-Fechner
[scores,allParams] = plotResults(logisticRLDirs{4},1,'fitStatsOnly',true,'Only120Trials',true);
figName = 'fig6B_logisticAbortRL_noWF_bestFit';
saveas(gcf,[FIGURE_DIR figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR figName '.png'],'-dpng','-r600')

%% Fig 7A model comparison (w/w.o. WF)
BICs = [];
nDataPoints = 160;

% Model comparison of w/w.o WF
% Load scores of best driftRL w WF
dirInd = 2;
[scores,allParams] = getScores(driftRLDirs{dirInd});
[~,minInd] = min(scores);
params = load([driftRLDirs{dirInd} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

% Load scores of best driftRL w.o WF
dirInd = 4;
[scores,allParams] = getScores(driftRLDirs{dirInd});
[~,minInd] = min(scores);
params = load([driftRLDirs{dirInd} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

% Load scores of best logisticAbortRL w WF
dirInd = 2;
[scores,allParams] = getScores(logisticRLDirs{dirInd});
[~,minInd] = min(scores);
params = load([logisticRLDirs{dirInd} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

% Load scores of best logisticAbortRL w WF
dirInd = 4;
[scores,allParams] = getScores(logisticRLDirs{dirInd});
[~,minInd] = min(scores);
params = load([logisticRLDirs{dirInd} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

xtick = [1:2 4:5];
xticklabels = {'DDM','DDM WF-','LOG(\tau)','LOG(\tau) WF-'};
figure;
elim_whitespace
bar(xtick,BICs)
set(gca,'xtick',xtick,'xticklabels',xticklabels,'fontsize',12,'fontweight','bold')
set(gcf,'Position',[10 10 1800 800])
ylabel('BIC','fontsize',20,'fontweight','bold')
figName = 'fig7A_BIC_with_without_WF';
saveas(gcf,[FIGURE_DIR '/' figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR '/' figName '.png'],'-dpng','-r600')

%% Fig 7B logistic vs. drift
BICs = [];
nDataPoints = 160;

% Load scores of best driftRL w WF
dirInd = 2;
[scores,allParams] = getScores(driftRLDirs{dirInd});
[~,minInd] = min(scores);
params = load([driftRLDirs{dirInd} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

% Load scores of best driftRL w WF
dirInd = 5;
[scores,allParams] = getScores(driftRLDirs{dirInd});
[~,minInd] = min(scores);
params = load([driftRLDirs{dirInd} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

% Load scores of best logisticAbortRL w WF
dirInd = 2;
[scores,allParams] = getScores(logisticRLDirs{dirInd});
[~,minInd] = min(scores);
params = load([logisticRLDirs{dirInd} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

% Load scores of best logisticAbortRL w WF
dirInd = 5;
[scores,allParams] = getScores(logisticRLDirs{dirInd});
[~,minInd] = min(scores);
params = load([logisticRLDirs{dirInd} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];

%{
% Load scores of best logisticAbortRL w WF
dirInd = 6;
[scores,allParams] = getScores(logisticRLDirs{dirInd});
[~,minInd] = min(scores);
params = load([logisticRLDirs{dirInd} '/' num2str(minInd) '/params.mat']); params=params.params;
curBIC = length(params)*log(nDataPoints) + 2*scores(minInd);
BICs = [BICs curBIC];
%}

%xtick = [1:2 4:6];
xtick = [1:2 4:5];
%xticklabels = {'DDM(\tau)','DDM(UCB)','LOG(\tau)','LOG(UCB)','LOG(\epsilon)'};
xticklabels = {'DDM(\tau)','DDM(UCB)','LOG(\tau)','LOG(UCB)'};
figure;
%elim_whitespace
bar(xtick,BICs - (min(BICs) - 10))
set(gca,'xtick',xtick,'xticklabels',xticklabels,'fontsize',20,'fontweight','bold')
ticks=get(gca,'YTick');%retrieve current ticks
ticks=ceil(ticks + min(BICs) - 10); % add back actual BIC values
%ticks=cellfun(@num2str,ticks,'UniformOutput',false);%convert to cellstr
set(gca,'YTickLabels',ticks)%set new tick labels
set(gcf,'Position',[10 10 1800 800])
ylabel('BIC','fontsize',20,'fontweight','bold')
figName = 'fig7B_BIC_DDM_v_LOG';
saveas(gcf,[FIGURE_DIR '/' figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR '/' figName '.png'],'-dpng','-r600')

%% Alternative Fig 7B (using boxplots of multiple BICs)
BICs = zeros(20,5);
data_dir = 'publication/modeling/data_files/';

[scores,allParams] = getScores([data_dir 'DDM_' driftSimDataDirs{2}]);
BICs(:,1) = scores';

[scores,allParams] = getScores([data_dir 'DDM_' driftSimDataDirs{5}]);
BICs(:,2) = scores';

[scores,allParams] = getScores([data_dir 'logistic_' logisticSimDataDirs{2}]);
BICs(:,3) = scores';

[scores,allParams] = getScores([data_dir 'logistic_' logisticSimDataDirs{5}]);
BICs(:,4) = scores';

[scores,allParams] = getScores([data_dir 'logistic_' logisticSimDataDirs{6}]);
BICs(:,5) = scores';

xticklabels = {'DDM(\tau)','DDM(UCB)','LOG(\tau)','LOG(UCB)','LOG(\epsilon)'};
figure;
notBoxPlot(BICs)
ylabel('BIC','FontSize',20,'FontWeight','bold')
set(gca,'xticklabels',xticklabels,'FontSize',20,'FontWeight','bold')
set(gcf,'Position',[10 10 1200 1000])

%% Fig 8 - matching law
sideChoices = dataTable.SideChosen;
relativeValues = dataTable.RelativeValue;
PRvalues = dataTable.PRvalue;
FRvalues = dataTable.FRvalue;
valRatios = PRvalues./FRvalues;
PRvalFrac = PRvalues./(PRvalues + FRvalues);
lowVal = 0;
binWidth = .01;
maxVal = 1;
bins_L = lowVal:binWidth:(maxVal-binWidth);
bins_U = (lowVal+binWidth):binWidth:maxVal;
mid_bins = mean([bins_L; bins_U],1);
nbins = length(mid_bins);
pPR = zeros(1,nbins);
ratioPRFR = zeros(1,nbins);
for i=1:nbins
    relValInds = find(PRvalFrac >= bins_L(i) & PRvalFrac < bins_U(i));
    if (isempty(relValInds))
        pPR(i) = nan;
    else
        pPR(i) = sum(sideChoices(relValInds))/length(relValInds);
    end
    ratioPRFR(i) = pPR(i)/(1-pPR(i));
    
end

fn = @(b,x) b(4) + b(3)./(1 + exp(-b(1)*(x - b(2))));
beta0 = [1 .5 .5 0];
NLM = fitnlm(PRvalFrac,sideChoices,fn,beta0);
NLM2 = fitnlm(mid_bins,pPR,fn,beta0);

x = lowVal:binWidth:maxVal;
x2 = lowVal:(maxVal-lowVal)/1000:maxVal;
figure;

scatter(mid_bins,pPR,100,'k.')
hold on; plot(x,x)
%plot(x2,NLM.predict(x2'))
plot(x2,NLM2.predict(x2'))
xlabel('EoR_{PR} / (EoR_{PR} + EoR_{FR})','FontSize',20,'FontWeight','bold')
ylabel('P(PR)','FontSize',20,'FontWeight','bold')
legend({'data','Matching law','logistic fit'},'location','northwest','FontSize',20,'FontWeight','bold')
set(gcf,'Position',[10 10 1000 800])
figName = 'fig8_matching_law';
saveas(gcf,[FIGURE_DIR '/' figName '.fig'],'fig')
saveas(gcf,[FIGURE_DIR figName '.eps'],'epsc')
print([FIGURE_DIR '/' figName '.png'],'-dpng','-r600')

