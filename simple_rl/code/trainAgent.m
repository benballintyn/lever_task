function [Q,params,sessionTypes,allActions,allRewards] = trainAgent(agentType,actionSelectionMethod,agentParams,varargin)
p = inputParser;
isValidAgentType = @(x) any(strcmp(x,{'bandit','Qlearner_small','Qlearner_big'}));
isValidActionSelectionMethod = @(x) any(strcmp(x,{'e_greedy','softmax'}));
isFunc = @(x) isa(x,'function_handle');
isPositiveNumber = @(x) isnumeric(x) && x > 0;
addRequired(p,'agentType',isValidAgentType)
addRequired(p,'actionSelectionMethod',isValidActionSelectionMethod);
addRequired(p,'agentParams',@isstruct)
addParameter(p,'utilityFunc1',@(x) x,isFunc)
addParameter(p,'utilityFunc2',@(x) x,isFunc)
addParameter(p,'valueAdjustment','none',@ischar)
addParameter(p,'rewardType','raw',@ischar)
addParameter(p,'nTrials',500,isPositiveNumber)
addParameter(p,'nStates',200,isPositiveNumber);
addParameter(p,'nDays',16,isPositiveNumber);
parse(p,agentType,actionSelectionMethod,agentParams,varargin{:})
% Create params stucture
params.agentType = p.Results.agentType;
params.actionSelectionMethod = p.Results.actionSelectionMethod;
params.agentParams = p.Results.agentParams;
params.utilityFunc1 = p.Results.utilityFunc1;
params.utilityFunc2 = p.Results.utilityFunc2;
params.valueAdjustment = p.Results.valueAdjustment;
params.rewardType = p.Results.rewardType;
params.nTrials = p.Results.nTrials;
params.nStates = p.Results.nStates;
params.nDays = p.Results.nDays;

% Check that all necessary parameters are contained in agentParams
if (~any(strcmp('alpha',fieldnames(agentParams))))
    error('agentParams does not contain a learning rate alpha')
end
alpha = p.Results.agentParams.alpha;
switch p.Results.agentType
    case 'bandit'
        Q = cell(1,4);
        for i=1:4
            Q{i} = rand(1,2);
        end
    case 'Qlearner_small'
        if (~isfield(agentParams,'gamma'))
            error('agentParams is missing a discount factor gamma')
        else
            gamma = p.Results.agentParams.gamma;
        end
        Q = cell(1,4);
        for i=1:4
            Q{i} = rand(1,2);
        end
    case 'Qlearner_big'
        if (~isfield(agentParams,'gamma'))
            error('agentParams is missing a discount factor gamma')
        else
            gamma = p.Results.agentParams.gamma;
        end
        Q = cell(1,4);
        for i=1:4
            Q{i} = rand(p.Results.nStates,2);
        end
    otherwise
        error('agentType not recognized')
end
disp(['agentType               : ' p.Results.agentType])
disp(['nStates                 : ' num2str(size(Q{1},1))])
switch p.Results.actionSelectionMethod
    case 'e_greedy'
        if (~isfield(agentParams,'epsilon'))
            error('actionSelectionMethod is e_greedy but no epsilon exists in agentParams')
        else
            epsilon = p.Results.agentParams.epsilon;
        end
    case 'softmax'
        if (~isfield(agentParams,'temp'))
            error('actionSelectionMethod is softmax but no temp exists in agentParams')
        else
            temp = p.Results.agentParams.temp;
        end
    otherwise
        error('actionSelectionMethod not recognized')
end
disp(['actionSelectionMethod   : ' p.Results.actionSelectionMethod])
switch p.Results.valueAdjustment
    case 'UCB'
        if (~isfield(agentParams,'c'))
            error('valueAdjustmment is UCB but no c exists in agentParams')
        else
            c = p.Results.agentParams.c;
        end
end
disp(['valueAdjustment         : ' p.Results.valueAdjustment])
disp(['rewardType              : ' p.Results.rewardType])
disp(['nDays                   : ' num2str(p.Results.nDays)])
disp(['nTrials                 : ' num2str(p.Results.nTrials)])
disp(p.Results.agentParams)

% Pull some variables out of structure for performance
utilityFunc1 = p.Results.utilityFunc1;
utilityFunc2 = p.Results.utilityFunc2;
rewardType = p.Results.rewardType;
valueAdjustment = p.Results.valueAdjustment;
agentType = p.Results.agentType;
nDays = p.Results.nDays;
nTrials = p.Results.nTrials;
nStates = p.Results.nStates;

% Do training
for d = 1:nDays
    state = 1;
    newState = 1;
    sessionType = ceil(rand*4);
    sessionTypes(d) = sessionType;
    switch sessionType
        case 1 % 2xFR6
            Pl = 2; Ps = 6; SR = 3; LR = 6;
        case 2 % 2xFR12
            Pl = 2; Ps = 12; SR = 3; LR = 6;
        case 3 % 5xFR6
            Pl = 2; Ps = 6; SR = 3; LR = 15;
        case 4 % 5xFR12
            Pl = 2; Ps = 12; SR = 3; LR = 15;
    end
    utilityFuncValsSR = utilityFunc1(ones(1,nTrials)*Ps);
    utilityFuncValsPR = utilityFunc1(1:(nTrials+1));
    for j=1:length(utilityFuncValsSR)
        denomsSR(j) = utilityFunc2(utilityFuncValsSR(j));
    end
    for j=1:length(utilityFuncValsPR)
        denomsPR(j) = utilityFunc2(utilityFuncValsPR(j));
    end
    rewards = zeros(1,p.Results.nTrials);
    actions = zeros(1,p.Results.nTrials);
    switch valueAdjustment
        case 'UCB'
            action_count = zeros(1,2);
    end
    for i=1:nTrials
        % Choose action based on actionSelectionMethod and optional
        % valueAdjustment
        switch actionSelectionMethod
            case 'e_greedy'
                if (rand < epsilon)
                    action = ceil(rand*2);
                else
                    switch valueAdjustment
                        case 'none'
                            [~,action] = max(Q{sessionType}(state,:));
                        case 'UCB'
                            adjustedValues = Q{sessionType}(state,:) + c*sqrt(log(i)./action_count);
                            [~,action] = max(adjustedValues);
                            action_count(action) = action_count(action) + 1;
                        otherwise
                            error(['valueAdjustment: ' valueAdjustment ' not recognized'])
                    end
                end
            case 'softmax'
                switch valueAdjustment
                    case 'none'
                        [softVals,action] = mySoftmax(Q{sessionType}(state,:),temp);
                    case 'UCB'
                        adjustedValues = Q{sessionType}(state,:) + c*sqrt(log(i)./action_count);
                        [softvals,action] = mySoftmax(adjustedValues,temp);
                        action_count(action) = action_count(action) + 1;
                    otherwise
                        error(['valueAdjustment: ' valueAdjustment ' not recognized'])
                end
        end
        actions(i) = action;
        
        if (action == 1) % if chose action 1 (SR side)
            switch rewardType
                case 'raw'
                    r = SR;
                case 'divisive'
                    r = SR/denomsSR(i); %SR/utilityFunc2(utilityFuncValsSR(t));
                otherwise
                    error('rewardType not recognize')
            end
        else % if chose action 2 (PR side)
            switch rewardType
                case 'raw'
                    r = LR;
                case 'divisive'
                    r = LR/denomsPR(Pl); %LR/utilityFunc2(utilityFuncValsPR(Pl));
                otherwise
                    error('rewardType not recognized')
            end
            Pl = Pl + 1; % increment the PR side press cost
            switch agentType
                case 'bandit'
                    newState = 1;
                case 'Qlearner_small'
                    newState = 1;
                case 'Qlearner_big'
                    newState = state + 1;
            end
        end
        rewards(i) = r;
        if (newState > nStates)
            newState = nStates;
        end
        switch agentType
            case 'bandit'
                Q{sessionType}(state,action) = Q{sessionType}(state,action) + alpha*(r - Q{sessionType}(state,action));
            case 'Qlearner_small'
                Q{sessionType}(state,action) = Q{sessionType}(state,action) + alpha*(r + gamma*max(Q{sessionType}(newState,:)) - Q{sessionType}(state,action));
            case 'Qlearner_big'
                Q{sessionType}(state,action) = Q{sessionType}(state,action) + alpha*(r + gamma*max(Q{sessionType}(newState,:)) - Q{sessionType}(state,action));
        end
        state = newState;
    end
    allActions{d} = actions;
    allRewards{d} = rewards;
    disp(['Done with day #' num2str(d)])
end
end

