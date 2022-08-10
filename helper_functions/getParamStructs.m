function [agentParams,abortModelParams,forgettingParams] = getParamStructs(modelType,params,actionSelectionMethod,utilityFuncs,forgettingType)
% Parse parameters
agentParams.alpha = params(1);
switch actionSelectionMethod
    case 'e_greedy'
        agentParams.epsilon = params(2);
    case 'softmax'
        agentParams.temp = params(2);
    case 'UCB'
        agentParams.c = params(2);
end
if (any(strcmp(utilityFuncs,'ansUtilityFunc')))
    agentParams.ans_sigma = params(3);
    if (any(strcmp(modelType,{'driftRL','driftRL_valueUpdate'})))
        abortModelParams.drift_rate = params(4);
        abortModelParams.noise_amplitude = params(5);
    elseif (strcmp(modelType,'logisticAbortRL'))
        abortModelParams.temp = params(4);
        abortModelParams.offset = params(5);
    else
        error('modelType not recognized')
    end
    nextParamInd = 6;
else
    if (any(strcmp(modelType,{'driftRL','driftRL_valueUpdate'})))
        abortModelParams.drift_rate = params(3);
        abortModelParams.noise_amplitude = params(4);
    elseif (strcmp(modelType,'logisticAbortRL'))
        abortModelParams.temp = params(3);
        abortModelParams.offset = params(4);
    else
        error('modelType not recognized')
    end
    nextParamInd = 5;
end
switch forgettingType
    case 'none'
        if (length(params) > 6)
            error('Without forgetting, there should be at most 6 parameters')
        end
    case 'decayToInitialValues'
        forgettingParams = params(nextParamInd);
        if (length(params) > nextParamInd)
            error('There are extra parameters ???')
        end
    case 'decayToFreeParameter'
        forgettingParams = params(nextParamInd:(nextParamInd+1));
        if (length(params) > (nextParamInd+1))
            error('There are extra parameters ???')
        end
    otherwise
        error('forgettingType not recognized')
end
end

