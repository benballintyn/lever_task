% Set paths and filenames
clear all;
BASE_DIR = '~/phd/lever_task/publication/modeling/';
DATA_DIR = [BASE_DIR 'data_files/'];

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
    driftSimDataDirs{6} = 'bandit_e_greedy_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues_noAbortANS';
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

%% Run simulations
driftDirs2run = [6];

nruns = 100;

% First for DDM models
agentType = 'bandit';
initializationMethod = 'mean_reward';
forgettingType = 'decayToInitialValues';
scoreType = 'logprob_joint';
modelType = 'driftRL';
Only120Trials = true;
fullANS = false;
noAbortANS = true;
for i=1:length(driftDirs2run)
    CUR_DIR = [DATA_DIR 'DDM_' driftSimDataDirs{driftDirs2run(i)}];
    [scores,allParams] = getScores(driftRLDirs{driftDirs2run(i)});
    [~,minInd] = min(scores);
    params = load([driftRLDirs{driftDirs2run(i)} '/' num2str(minInd) '/params.mat']); params=params.params;
    
    if (exist(CUR_DIR,'dir'))
        s = getScores(CUR_DIR);
        curRunNum = length(s);
    else
        curRunNum = 0;
    end
    % determine actionSelectionMethod from filename
    if (contains(driftSimDataDirs{driftDirs2run(i)},'softmax'))
        actionSelectionMethod = 'softmax';
    elseif (contains(driftSimDataDirs{driftDirs2run(i)},'UCB'))
        actionSelectionMethod = 'UCB';
    elseif (contains(driftSimDataDirs{driftDirs2run(i)},'e_greedy'))
        actionSelectionMethod = 'e_greedy';
    else
        error('No action selection method identified')
    end
    
    % determine utility functions from filename
    if (contains(driftSimDataDirs{driftDirs2run(i)},'ansUtilityFunc'))
        utilityFunc1 = 'ansUtilityFunc';
        if (contains(driftSimDataDirs{driftDirs2run(i)},'pressUtilityFunc'))
            utilityFunc2 = 'pressUtilityFunc';
        else
            utilityFunc2 = '';
        end
    elseif (contains(driftSimDataDirs{driftDirs2run(i)},'pressUtilityFunc'))
        utilityFunc1 = '';
        utilityFunc2 = 'pressUtilityFunc';
    else
        utilityFunc1 = '';
        utilityFunc2 = '';
    end
    
    parfor j=1:nruns
        [score] = run_param_set(params,CUR_DIR,agentType,...
        actionSelectionMethod,initializationMethod,utilityFunc1,utilityFunc2,forgettingType,...
        scoreType,modelType,curRunNum+j,'Only120Trials',Only120Trials,'fullANS',fullANS,'noAbortANS',noAbortANS);
    end
end

%% Next for logistic models
logisticDirs2run = [2 5 6];
nruns = 20;

agentType = 'bandit';
initializationMethod = 'mean_reward';
forgettingType = 'decayToInitialValues';
scoreType = 'logprob_joint';
modelType = 'logisticAbortRL';
Only120Trials = true;
fullANS = false;
noAbortANS = true;
for i=1:length(logisticDirs2run)
    CUR_DIR = [DATA_DIR 'logistic_' logisticSimDataDirs{logisticDirs2run(i)}];
    [scores,allParams] = getScores(logisticRLDirs{logisticDirs2run(i)});
    [~,minInd] = min(scores);
    params = load([logisticRLDirs{logisticDirs2run(i)} '/' num2str(minInd) '/params.mat']); params=params.params;
    
    if (exist(CUR_DIR,'dir'))
        s = getScores(CUR_DIR);
        curRunNum = length(s);
    else
        curRunNum = 0;
    end
    % determine actionSelectionMethod from filename
    if (contains(logisticSimDataDirs{logisticDirs2run(i)},'softmax'))
        actionSelectionMethod = 'softmax';
    elseif (contains(logisticSimDataDirs{logisticDirs2run(i)},'UCB'))
        actionSelectionMethod = 'UCB';
    elseif (contains(logisticSimDataDirs{logisticDirs2run(i)},'e_greedy'))
        actionSelectionMethod = 'e_greedy';
    else
        error('No action selection method identified')
    end
    
    % determine utility functions from filename
    if (contains(logisticSimDataDirs{logisticDirs2run(i)},'ansUtilityFunc'))
        utilityFunc1 = 'ansUtilityFunc';
        if (contains(logisticSimDataDirs{logisticDirs2run(i)},'pressUtilityFunc'))
            utilityFunc2 = 'pressUtilityFunc';
        else
            utilityFunc2 = '';
        end
    elseif (contains(logisticSimDataDirs{logisticDirs2run(i)},'pressUtilityFunc'))
        utilityFunc1 = '';
        utilityFunc2 = 'pressUtilityFunc';
    else
        utilityFunc1 = '';
        utilityFunc2 = '';
    end
    
    parfor j=1:nruns
        [score] = run_param_set(params,CUR_DIR,agentType,...
        actionSelectionMethod,initializationMethod,utilityFunc1,utilityFunc2,forgettingType,...
        scoreType,modelType,curRunNum+j,'Only120Trials',Only120Trials,'fullANS',fullANS,'noAbortANS',noAbortANS);
    end
end

%%
