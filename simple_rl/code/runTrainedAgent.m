function [states,actions,rewards,numLR] = runTrainedAgent(Q,agentType,actionSelectionMethod,...
                                                        agentParams,nTrials,nRepetitions,varargin)
p = inputParser;
isValidAgentType = @(x) any(strcmp(x,{'bandit','Qlearner_small','Qlearner_big'}));
isValidActionSelectionMethod = @(x) any(strcmp(x,{'e_greedy','softmax'}));
isValidValueAdjustment = @(x) any(strcmp(x,{'none','UCB'}));
isValidSessionType = @(x) any(strcmp(x,{'2xFR6','2xFR12','5xFR6','5xFR12'}));
isValidRewardType = @(x) any(strcmp(x,{'raw','divisive'}));
isPositiveNumber = @(x) isnumeric(x) && x > 0;
addRequired(p,'Q',@isnumeric)
addRequired(p,'agentType',isValidAgentType)
addRequired(p,'actionSelectionMethod',isValidActionSelectionMethod)
addRequired(p,'agentParams',@isstruct)
addRequired(p,'nTrials',isPositiveNumber)
addRequired(p,'nRepetitions',isPositiveNumber)
addRequired(p,'sessionType',isValidSessionType)
addParameter(p,'useANS',false,@islogical)
addParameter(p,'valueAdjustment','none',isValidValueAdjustment)
addParameter(p,'rewardType','raw',isValidRewardType)
parse(p,Q,agentType,actionSelectionMethod,agentParams,nTrials,nRepetitions,varargin{:})
switch p.Results.agentType
    case 'bandit'
        alpha = p.Results.agentParams.alpha;
    case 'Qlearner_small'
        alpha = p.Results.agentParams.alpha;
        gamma = p.Results.agentParams.gamma;
    case 'Qlearner_big'
        alpha = p.Results.agentParams.alpha;
        gamma = p.Results.agentParams.gamma;
end
switch p.Results.actionSelectionMethod
    case 'e_greedy'
        epsilon = p.Results.agentParams.epsilon;
    case 'softmax'
        temp = p.Results.agentParams.temp;
end
switch p.Results.valueAdjustment
    case 'UCB'
        c = p.Results.agentParams.c;
end
if (p.Results.useANS)
    if (~isfield(p.Results.agentParams,'ans_sigma'))
        error('ans_sigma is not a field of agentParams')
    else
        ans_sigma = p.Results.agentParams.ans_sigma;
        ansUtilityFunc = @(x) lognrnd(log(x) - (ans_sigma^2)/2,ans_sigma);
    end
else
    ansUtilityFunc = @(x) x;
end
% pull variables out for better performance
nRepetitions = p.Results.nRepetitions;
nTrials = p.Results.nTrials;
actionSelectionMethod = p.Results.actionSelectionMethod;
valueAdjustment = p.Results.valueAdjustment;
rewardType = p.Results.rewardType;
sessionType = p.Results.sessionType;

% Initialize outputs
states = zeros(nRepetitions,nTrials);
actions = zeros(nRepetitions,nTrials);
rewards = zeros(nRepetitions,nTrials);
numLR = zeros(1,nRepetitions);
% Run agent through task nRepetitions times
for i=1:nRepetitions
    switch sessionType
        case '2xFR6'
            Pl = 2; Ps = 6; SR = 3; LR = 6;
        case '2xFR12'
            Pl = 2; Ps = 12; SR = 3; LR = 6;
        case '5xFR6'
            Pl = 2; Ps = 6; SR = 3; LR = 15;
        case '5xFR12'
            Pl = 2; Ps = 12; SR = 3; LR = 15;
    end
    switch valueAdjustment
        case 'UCB'
            action_count = zeros(1,2);
    end
    state = 1;
    ANS_valsSR = ansUtilityFunc(ones(1,nTrials)*Ps);
    ANS_valsPR = ansUtilityFunc(1:(nTrials+1));
    for j=1:nTrials
        % Choose action based on actionSeletionMethod and valueAdjustment
        switch actionSelectionMethod
            case 'e_greedy'
                if (rand < epsilon)
                    action = ceil(rand*2);
                else
                    switch valueAdjustment
                        case 'none'
                            [~,action] = max(Q(state,:));
                        case 'UCB'
                            adjustedValues = Q(state,:) + c*sqrt(log(j)./action_count);
                            [~,action] = max(adjustedValues);
                            action_count(action) = action_count(action) + 1;
                        otherwise
                            error(['valueAdjustment: ' valueAdjustment ' not recognized'])
                    end
                end
            case 'softmax'
                switch valueAdjustment
                    case 'none'
                        [softVals,action] = mySoftmax(Q(state,:),temp);
                    case 'UCB'
                        adjustedValues = Q(state,:) + c*sqrt(log(i)./action_count);
                        [softvals,action] = mySoftmax(adjustedValues,temp);
                        action_count(action) = action_count(action) + 1;
                    otherwise
                        error(['valueAdjustment: ' valueAdjustment ' not recognized'])
                end
        end
        actions(i,j) = action;
        
        % Determine reward based on rewardType and ANS
        if (action == 1) % if chose action 1 (SR side)
            switch rewardType
                case 'raw'
                    r = SR;
                case 'divisive'
                    r = SR/ANS_valsSR(i); %SR/utilityFunc2(utilityFuncValsSR(t));
                otherwise
                    error('rewardType not recognize')
            end
        else % if chose action 2 (PR side)
            switch rewardType
                case 'raw'
                    r = LR;
                case 'divisive'
                    r = LR/ANS_valsPR(Pl); %LR/utilityFunc2(utilityFuncValsPR(Pl));
                otherwise
                    error('rewardType not recognized')
            end
            Pl = Pl + 1; % increment the PR side press cost
            switch agentType
                case 'Qlearner_big'
                    state = state + 1;
            end
        end
        rewards(i,j) = r;
        states(i,j) = state;
    end
    numLR(i) = Pl - 1;
end
end

