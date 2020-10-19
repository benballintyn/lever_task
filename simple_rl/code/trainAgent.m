function [Q,params] = trainAgent(agentType,actionSelectionMethod,agentParams,varargin)
p = inputParser;
isValidAgentType = @(x) any(strcmp(x,{'bandit','Qlearner_small','Qlearner_big'}));
isValidActionSelectionMethod = @(x) any(strcmp(x,{'e_greedy','softmax'}));
isFunc = @(x) isa(x,'function_handle');
addRequired(p,'agentType',isValidAgentType)
addRequired(p,'actionSelectionMethod',isValidActionSelectionMethod);
addRequired(p,'agentParams',@isstruct)
addParameter(p,'utilityFunc1',@(x) x,isFunc)
addParameter(p,'utilityFunc2',@(x) x,isFunc)
addParameter(p,'valueAdjustment','none',@ischar)
addParameter(p,'rewardType','raw',@ischar)
addParameter(p,'nTrials',500,@isnumeric)
addParameter(p,'nStates',200,@isnumeric);
addParameter(p,'nDays',16,@isnumeric);
parse(p,agentType,actionSelectionMethod,agentParams,varargin{:})

% Check that all necessary parameters are contained in agentParams
if (~any(strcmp('alpha',fieldnames(agentParams))))
    error('agentParams does not contain a learning rate alpha')
end
switch p.Results.agentType
    case 'bandit'
        Q = cell(1,4);
        for i=1:4
            Q{i} = rand(1,2);
        end
    case 'Qlearner_small'
        if (~isfield(agentParams,'gamma'))
            error('agentParams is missing a discount factor gamma')
        end
        Q = cell(1,4);
        for i=1:4
            Q{i} = rand(1,2);
        end
    case 'Qlearner_big'
        if (~isfield(agentParams,'gamma'))
            error('agentParams is missing a discount factor gamma')
        end
        Q = cell(1,4);
        for i=1:4
            Q{i} = rand(p.Results.nStates,2);
        end
end
switch p.Results.actionSelectionMethod
    case 'e_greedy'
        if (~isfield(agentParams,'epsilon'))
            error('actionSelectionMethod is e_greedy but no epsilon exists in agentParams')
        end
    case 'softmax'
        if (~isfield(agentParams,'temp'))
            error('actionSelectionMethod is softmax but no temp exists in agentParams')
        end
end
switch p.Results.valueAdjustment
    case 'UCB'
        if (~isfield(agentParams,'c'))
            error('valueAdjustmment is UCB but no c exists in agentParams')
        else
            action_count = zeros(1,2);
        end
end
utilityFunc1 = p.Results.utilityFunc1;
utilityFunc2 = p.Results.utilityFunc2;
for d = 1:p.Results.nDays
    state = 1;
    sessionType = ceil(rand*4);
    switch sessionType
        case 1 % 2xFR6
            Pl = 2; Ps = 6; SR = 3; PR = 6;
        case 2 % 2xFR12
            Pl = 2; Ps = 12; SR = 3; PR = 6;
        case 3 % 5xFR6
            Pl = 2; Ps = 6; SR = 3; PR = 15;
        case 4 % 5xFR12
            Pl = 2; Ps = 12; SR = 3; PR = 15;
    end
    utilityFuncValsSR = utilityFunc1(ones(1,nTrials)*Ps);
    utilityFuncValsPR = utilityFunc1(1:(nTrials+1));
    for j=1:length(utilityFuncValsSR)
        denomsSR(j) = utilityFunc2(utilityFuncValsSR(j));
    end
    for j=1:length(utilityFuncValsPR)
        denomsPR(j) = utilityFunc2(utilityFuncValsPR(j));
    end
    nPresses = 0;
    rewards = zeros(1,p.Results.nTrials);
    actions = zeros(1,p.Results.nTrials);
    presses = zeros(1,p.Results.nTrials);
    for i=1:p.Results.nTrials
        % Choose action based on actionSelectionMethod and optional
        % valueAdjustment
        switch actionSelectionMethod
            case 'e_greedy'
                if (rand < agentParams.epsilon)
                    action = ceil(rand*2);
                else
                    switch p.Results.valueAdjustment
                        case 'none'
                            [~,action] = max(Q(state,:));
                        case 'UCB'
                            adjustedValues = Q(state,:) + agentParams.c*sqrt(log(i)./action_count);
                            [~,action] = max(adjustedValues);
                        otherwise
                            error(['valueAdjustment: ' p.Results.valueAdjustment ' not recognized'])
                    end
                end
            case 'softmax'
                switch p.Results.valueAdjustment
                    case 'none'
                        [softVals,action] = mySoftmax(Q(state,:),temp);
                    case 'UCB'
                        adjustedValues = Q(state,:) + agentParams.c*sqrt(log(i)./action_count);
                        [softvals,action] = mySoftmax(adjustedValues,temp);
                    otherwise
                        error(['valueAdjustment: ' p.Results.valueAdjustment ' not recognized'])
                end
        end
        action_count(action) = action_count(action) + 1;
        actions(t) = action;
        
        if (action == 1) % if chose action 1 (SR side)
            r = SR/denomsSR(t); %SR/utilityFunc2(utilityFuncValsSR(t));
            rewards(t) = SR;
            nPresses = nPresses + Ps;
            nS(t) = nS(t) + 1;
            presses(t) = Ps;
        else % if chose action 2 (PR side)
            curLRtrial = curLRtrial+1;
            r = LR/denomsPR(t); %LR/utilityFunc2(utilityFuncValsPR(Pl));
            rewards(t) = LR;
            nPresses = nPresses + Pl;
            presses(t) = Pl;
            Pl = Pl + 1; % increment the PR side press cost
            nL(t) = nL(t) + 1;
        end
    end
end
end

