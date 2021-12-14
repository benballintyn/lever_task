agentType = 'bandit';
actionSelectionMethod = 'softmax';
utilityFunc1 = 'ansUtilityFunc';
utilityFunc2 = '';
initializationMethod = 'mean_reward';
forgettingType = 'decayToInitialValues';
scoreType = 'logprob_joint';
modelType = 'logisticAbortRL';
Only120Trials = true;
fullANS = true;

N_RUNS = 2000;

[status,SYSTEM_NAME] = system('hostname');
if (~status)
    switch strtrim(SYSTEM_NAME)
        case 'silmaril'
            savedir = ['/media/ben/Varda/phd/lever_task/' modelType '/results/Only120Trials/' scoreType];
        case 'miller-lab-ubuntu2'
            savedir = ['/media/ben/Manwe/phd/lever_task/' modelType '/results/Only120Trials/' scoreType];
        otherwise
            error('System hostname not recognized')
    end
end

[xmin] = lever_task_surrogate_optim(agentType,actionSelectionMethod,...
utilityFunc1,utilityFunc2,initializationMethod,forgettingType,scoreType,modelType,...
savedir,N_RUNS,'Only120Trials',Only120Trials,'fullANS',fullANS);
