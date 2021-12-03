% runGridSearch
agentType = 'bandit';
actionSelectionMethod = 'softmax';
utilityFunc1 = '';
utilityFunc2 = '';
initializationMethod = 'mean_reward';
forgettingType = 'decayToInitialValues';
scoreType = 'logprob_joint';
modelType = 'driftRL';
useOnly120Trials = true;
driftType = 'value_based_drift';

savedir = ['/media/ben/Manwe/phd/lever_task/' modelType '/results/Only120Trials/' scoreType];
if (~exist(savedir,'dir'))
    mkdir(savedir)
end

[scores] = grid_search(agentType,actionSelectionMethod,...
    utilityFunc1,utilityFunc2,initializationMethod,forgettingType,scoreType,...
    modelType,savedir,'Only120Trials',useOnly120Trials,'driftType',driftType);