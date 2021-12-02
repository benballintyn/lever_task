agentType = 'bandit';
actionSelectionMethod = 'softmax';
utilityFunc1 = 'ansUtilityFunc';
utilityFunc2 = 'pressUtilityFunc';
initializationMethod = 'mean_reward';
forgettingType = 'decayToInitialValues';
scoreType = 'logprob_joint';
modelType = 'driftRL';
savedir = ['/media/ben/Varda/phd/lever_task/' modelType '/results/Only120Trials/' scoreType];

[xmin] = lever_task_surrogate_optim(agentType,actionSelectionMethod,...
utilityFunc1,utilityFunc2,initializationMethod,forgettingType,scoreType,modelType,...
savedir,2000,'Only120Trials',true);
