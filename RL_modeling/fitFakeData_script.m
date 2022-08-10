agentType = 'bandit';
actionSelectionMethod = 'softmax';
utilityFunc1 = 'ansUtilityFunc';
utilityFunc2 = '';
initializationMethod = 'mean_reward';
forgettingType = 'decayToInitialValues';
modelType = 'logisticAbortRL';
maxFunEvals = 1000;
fullANS = false;
noAbortANS = true;
fakeData_type = 'bandit_softmax_ansUtilityFunc_initialization_mean_reward_forgettingType_decayToInitialValues_noAbortANS';
modelSpecStr = genModelSpecStr(agentType,actionSelectionMethod,utilityFunc1,utilityFunc2,initializationMethod,forgettingType,fullANS,noAbortANS);
for i=1:100
    basedir = ['/media/ben/Manwe/phd/lever_task/fitFakeData/logisticAbortRL/' fakeData_type '/' num2str(i) '/' modelType];
    if (strcmp(modelType,'driftRL'))
        basedir = [basedir '/value_based_drift'];
    end
    fakeDataDir = ['/media/ben/Manwe/phd/lever_task/fake_data/logisticAbortRL/' fakeData_type '/' num2str(i)];
    [xmin] = fitFakeData(agentType,actionSelectionMethod,...
        utilityFunc1,utilityFunc2,initializationMethod,forgettingType,...
        modelType,basedir,maxFunEvals,fakeDataDir,'fullANS',fullANS,'noAbortANS',noAbortANS);
end

modelType = 'driftRL';
for i=1:100
    basedir = ['/media/ben/Manwe/phd/lever_task/fitFakeData/logisticAbortRL/' fakeData_type '/' num2str(i) '/' modelType];
    if (strcmp(modelType,'driftRL'))
        basedir = [basedir '/value_based_drift'];
    end
    fakeDataDir = ['/media/ben/Manwe/phd/lever_task/fake_data/logisticAbortRL/' fakeData_type '/' num2str(i)];
    [xmin] = fitFakeData(agentType,actionSelectionMethod,...
        utilityFunc1,utilityFunc2,initializationMethod,forgettingType,...
        modelType,basedir,maxFunEvals,fakeDataDir,'fullANS',fullANS,'noAbortANS',noAbortANS);
end