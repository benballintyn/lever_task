% runGridSearch
agentType = 'bandit';
actionSelectionMethod = 'UCB';
utilityFunc1 = 'ansUtilityFunc';
utilityFunc2 = '';
initializationMethod = 'mean_reward';
forgettingType = 'decayToInitialValues';
scoreType = 'logprob_joint';
modelType = 'logisticAbortRL';
useOnly120Trials = true;
driftType = 'value_based_drift';
useFullANS = false;
noAbortANS = true;

savedir = ['/media/ben/Manwe/phd/lever_task/' modelType '/results/Only120Trials/' scoreType];
if (~exist(savedir,'dir'))
    mkdir(savedir)
end

[scores] = grid_search(agentType,actionSelectionMethod,...
    utilityFunc1,utilityFunc2,initializationMethod,forgettingType,scoreType,...
    modelType,savedir,'Only120Trials',useOnly120Trials,'driftType',driftType,'fullANS',useFullANS,'noAbortANS',noAbortANS);
