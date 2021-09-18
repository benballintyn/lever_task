function [lb,ub] = getParamBounds(agentType,actionSelectionMethod,utilityFuncs,forgettingType,modelType)
isValidAgentType = @(x) ismember(x,{'bandit','Qlearner'});
isValidActionSelectionMethod = @(x) ismember(x,{'e_greedy','softmax'});
areValidUtilityFuncs = @(x) isempty(x) || (length(x) <= 2)
isValidForgettingType = @(x) ismember(x,{'none','decayToInitialValues','decayToFreeParameter'})
isValidModelType = @(x) ismember(x,{'driftRL','driftRL_valueUpdate','logisticAbortRL'})
p = inputParser;
addRequired(p,'agentType',isValidAgentType)
addRequired(p,'actionSelectionMethod',isValidActionSelectionMethod)
addRequired(p,'utilityFuncs',areValidUtilityFuncs)
addRequired(p,'forgettingType',isValidForgettingType)
addRequired(p,'modelType',isValidModelType)
parse(p,agentType,actionSelectionMethod,utilityFuncs,forgettingType)

switch p.Results.agentType
    case 'bandit'
        lb(1) = 0; ub(1) = 1; % alpha - learning rate
        switch p.Results.actionSelectionMethod
            case 'e_greedy'
                lb(2) = 1e-4; ub(2) = 1; % epsilon - exploration rate
            case 'softmax'
                lb(2) = 1e-12; ub(2) = 1; % temp - exploration rate
        end
        nextParamInd = 3;
        if (any(strcmp(p.Results.utilityFuncs,'ansUtilityFunc')))
            lb(3) = 0; ub(3) = 2; % ANS sigma
            nextParamInd = 4;
        end
        if (any(strcmp(p.Results.modelType,{'driftRL','driftRL_valueUpdate'})))
            lb(nextParamInd) = -100; ub(nextParamInd) = 100; % drift rate coeff
            lb(nextParamInd+1) = 0; ub(nextParamInd+1) = 50; % drift noise scale
            nextParamInd = nextParamInd + 2;
        elseif (strcmp(p.Results.modelType,'logisticAbortRL'))
            lb(nextParamInd) = 1e-12; ub(nextParamInd) = 10; % logisticAbort process temp
            lb(nextParamInd+1) = -5; ub(nextParamInd+1) = 5; % logisticAbort process offset
            nextParamInd = nextParamInd + 2;
        else
            error('modelType not recognized')
        end
    case 'Qlearner'
        lb(1:2) = 0; ub(1:2) = 1; % alpha/gamma - learning rate / discount
        switch p.Results.actionSelectionMethod
            case 'e_greedy'
                lb(3) = 1e-4; ub(3) = 1; % epsilon - exploration rate
            case 'softmax'
                lb(3) = 1e-12; ub(3) = 1; % temp - exploration rate
        end
        nextParamInd = 4;
        if (any(strcmp(utilityFuncs,'ansUtilityFunc')))
            lb(4) = 0; ub(4) = 2; % ANS sigma
            nextParamInd = 5;
        end

        if (any(strcmp(p.Results.modelType,{'driftRL','driftRL_valueUpdate'})))
            lb(nextParamInd) = -100; ub(nextParamInd) = 100; % drift rate coeff
            lb(nextParamInd+1) = 0; ub(nextParamInd+1) = 50; % drift noise scale
            nextParamInd = nextParamInd + 2;
        elseif (strcmp(p.Results.modelType,'logisticAbortRL'))
            lb(nextParamInd) = 1e-12; ub(nextParamInd) = 10; % logisticAbort process temp
            lb(nextParamInd+1) = -5; ub(nextParamInd+1) = 5; % logisticAbort process offset
            nextParamInd = nextParamInd + 2;
        else
            error('modelType not recognized')
        end
end

switch p.Results.forgettingType
    case 'none'
    case 'decayToInitialValues'
        lb(nextParamInd) = 0; ub(nextParamInd) = 1; % alphaF (forgetting process 'learning rate')
    case 'decayToFreeParameter'
        lb(nextParamInd) = 0; ub(nextParamInd) = 1; % alphaF (forgetting process 'learning rate')
        lb(nextParamInd+1) = -2; ub(nextParamInd) = 2; %forgettingTargetCoeff
    otherwise
        error('forgettingType not recognized')
end
end
