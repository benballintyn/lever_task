function [lb,ub,paramNames] = getParamBounds(agentType,actionSelectionMethod,utilityFuncs,forgettingType,modelType)
isValidAgentType = @(x) ismember(x,{'bandit','Qlearner'});
isValidActionSelectionMethod = @(x) ismember(x,{'e_greedy','softmax'});
areValidUtilityFuncs = @(x) isempty(x) || (length(x) <= 2);
isValidForgettingType = @(x) ismember(x,{'none','decayToInitialValues','decayToFreeParameter'});
isValidModelType = @(x) ismember(x,{'driftRL','driftRL_valueUpdate','logisticAbortRL'});
p = inputParser;
addRequired(p,'agentType',isValidAgentType)
addRequired(p,'actionSelectionMethod',isValidActionSelectionMethod)
addRequired(p,'utilityFuncs',areValidUtilityFuncs)
addRequired(p,'forgettingType',isValidForgettingType)
addRequired(p,'modelType',isValidModelType)
parse(p,agentType,actionSelectionMethod,utilityFuncs,forgettingType,modelType)

paramNames = {};
switch p.Results.agentType
    case 'bandit'
        lb(1) = 0; ub(1) = 1; % alpha - learning rate
        paramNames{1} = 'alpha';
        switch p.Results.actionSelectionMethod
            case 'e_greedy'
                lb(2) = 1e-4; ub(2) = 1; % epsilon - exploration rate
                paramNames{2} = 'epsilon';
            case 'softmax'
                lb(2) = 1e-12; ub(2) = 1; % temp - exploration rate
                paramNames{2} = 'temperature';
            case 'UCB'
                lb(2) = 1e-12; ub(2) = 10; % c - exploration magnitude
                paramNames{2} = 'c';
            otherwise
                error(['Action selection method (' num2str(p.Results.actionSelectionMethod) ') not recognized'])
        end
        nextParamInd = 3;
        if (any(strcmp(p.Results.utilityFuncs,'ansUtilityFunc')))
            lb(3) = 0; ub(3) = 2; % ANS sigma
            paramNames{3} = 'ANS sigma';
            nextParamInd = 4;
        end
        if (any(strcmp(p.Results.modelType,{'driftRL','driftRL_valueUpdate'})))
            lb(nextParamInd) = -100; ub(nextParamInd) = 100; % drift rate coeff
            paramNames{nextParamInd} = 'drift rate coefficient';
            lb(nextParamInd+1) = 0; ub(nextParamInd+1) = 50; % drift noise scale
            paramNames{nextParamInd+1} = 'drift noise scale';
            nextParamInd = nextParamInd + 2;
        elseif (strcmp(p.Results.modelType,'logisticAbortRL'))
            lb(nextParamInd) = 1e-12; ub(nextParamInd) = 10; % logisticAbort process temp
            paramNames{nextParamInd} = 'logisticAbort temperature';
            lb(nextParamInd+1) = -5; ub(nextParamInd+1) = 5; % logisticAbort process offset
            paramNames{nextParamInd+1} = 'logisticAbort offset';
            nextParamInd = nextParamInd + 2;
        else
            error('modelType not recognized')
        end
    case 'Qlearner'
        lb(1:2) = 0; ub(1:2) = 1; % alpha/gamma - learning rate / discount
        paramNames{1} = 'alpha';
        paramNames{2} = 'gamma';
        switch p.Results.actionSelectionMethod
            case 'e_greedy'
                lb(3) = 1e-4; ub(3) = 1; % epsilon - exploration rate
                paramNames{3} = 'epsilon';
            case 'softmax'
                lb(3) = 1e-12; ub(3) = 1; % temp - exploration rate
                paramNames{3} = 'temperature';
        end
        nextParamInd = 4;
        if (any(strcmp(utilityFuncs,'ansUtilityFunc')))
            lb(4) = 0; ub(4) = 2; % ANS sigma
            paramNames{4} = 'ANS sigma';
            nextParamInd = 5;
        end

        if (any(strcmp(p.Results.modelType,{'driftRL','driftRL_valueUpdate'})))
            lb(nextParamInd) = -100; ub(nextParamInd) = 100; % drift rate coeff
            paramNames{nextParamInd} = 'drift rate coefficient';
            lb(nextParamInd+1) = 0; ub(nextParamInd+1) = 50; % drift noise scale
            paramNames{nextParamInd + 1} = 'drift noise scale';
            nextParamInd = nextParamInd + 2;
        elseif (strcmp(p.Results.modelType,'logisticAbortRL'))
            lb(nextParamInd) = 1e-12; ub(nextParamInd) = 10; % logisticAbort process temp
            paramNames{nextParamInd} = 'logisticAbort temperature';
            lb(nextParamInd+1) = -5; ub(nextParamInd+1) = 5; % logisticAbort process offset
            paramNames{nextParamInd + 1} = 'logisticAbort offset';
            nextParamInd = nextParamInd + 2;
        else
            error('modelType not recognized')
        end
end

switch p.Results.forgettingType
    case 'none'
    case 'decayToInitialValues'
        lb(nextParamInd) = 0; ub(nextParamInd) = 1; % alphaF (forgetting process 'learning rate')
        paramNames{nextParamInd} = 'alphaF';
    case 'decayToFreeParameter'
        lb(nextParamInd) = 0; ub(nextParamInd) = 1; % alphaF (forgetting process 'learning rate')
        paramNames{nextParamInd} = 'alphaF';
        lb(nextParamInd+1) = -2; ub(nextParamInd) = 2; %forgettingTargetCoeff
        paramNames{nextParamInd} = 'forgettingTargetCoeff';
    otherwise
        error(['forgettingType (' p.Results.forgettingType ') not recognized'])
end
end
