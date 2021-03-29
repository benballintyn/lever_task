function [Q,actions,rewards,trialEoRs,curLRtrial] = testSimpleAgent(sessionType,nTrials,agentType,actionSelectionMethod,agentParams,varargin)
% set up inputParser object and parse inputs
positiveNoInfCheck = @(x) x > 0 && ~isinf(x);
isaFunc = @(x) isa(x,'function_handle');
p = inputParser;
addRequired(p,'sessionType',@ischar) % 2xFR6, 2xFR12, 5xFR6, 5xFR12
addRequired(p,'nTrials',positiveNoInfCheck) %  # of trials to run per agent
addRequired(p,'agentType',@ischar) % currently: 'bandit' or 'Qlearner'
addRequired(p,'actionSelectionMethod',@ischar) % currently: 'e_greedy' or 'softmax'
addRequired(p,'agentParams',@isstruct) % contains the parameters of the agent (alpha,gamma,epsilon/temp)
addParameter(p,'utilityFunc1',@(x) x,isaFunc) % optional utility function to pass lever presses through
addParameter(p,'utilityFunc2',@(x) x,isaFunc) % optional utility function to pass the results of utilityFunc1 through
addParameter(p,'initialization','random',@ischar)
parse(p,sessionType,nTrials,agentType,actionSelectionMethod,agentParams,varargin{:})

nTrials = p.Results.nTrials;
sessionType = p.Results.sessionType;
agentType = p.Results.agentType;
actionSelectionMethod = p.Results.actionSelectionMethod;
agentParams = p.Results.agentParams;
utilityFunc1 = p.Results.utilityFunc1;
utilityFunc2 = p.Results.utilityFunc2;

switch agentType
    case 'bandit'
        alpha = p.Results.agentParams.alpha;
        %disp(['alpha = ' num2str(alpha)])
    case 'Qlearner'
        alpha = agentParams.alpha;
        gamma = agentParams.gamma;
        %disp(['alpha = ' num2str(alpha) '  gamma = ' num2str(gamma)])
    otherwise
        error('agentType not recognized')
end
switch actionSelectionMethod
    case 'e_greedy'
        epsilon = p.Results.agentParams.epsilon;
        %disp(['epsilon = ' num2str(epsilon)])
    case 'softmax'
        temp = p.Results.agentParams.temp;
       % disp(['temp = ' num2str(temp)])
    otherwise
        error('actionSelectionMethod not recognized')
end
switch sessionType
    case '2xFR6' % 2xFR6
        SR = 3; LR = 6; Ps = 6; sessInd = 1;
    case '2xFR12' % 2xFR12
        SR = 3; LR = 6; Ps = 12; sessInd = 2;
    case '5xFR6' % 5xFR6
        SR = 3; LR = 15; Ps = 6; sessInd = 3;
    case '5xFR12' % 5xFR12
        SR = 3; LR = 15; Ps = 12; sessInd = 4;
end
Pl = 2; % initialize PR side lever press cost
switch p.Results.initialization
    case 'random'
        Q = zeros(nTrials,2); Q(1,:) = rand(1,2); % initialize a 2-element Q table with uniform random values (0,1)
    case 'trained'
        Q = zeros(nTrials,2); Q(1,1) = SR/Ps; Q(1,2) = LR/Pl;
    case 'mean_reward'
        Q = zeros(nTrials,2); Q(1,1) = mean([SR/Ps LR/Pl]); Q(1,2) = mean([SR/Ps LR/Pl]);
end
nPresses = 0;
rewards = zeros(1,nTrials);
actions = zeros(1,nTrials);
presses = zeros(1,nTrials);
curLRtrial = 0;
utilityFuncValsSR = utilityFunc1(ones(1,nTrials)*Ps);
utilityFuncValsPR = utilityFunc1(1:(nTrials+1));
for t=2:nTrials
    % Choose action based on actionSelectionMethod
    switch actionSelectionMethod
        case 'e_greedy'
            if (rand < epsilon)
                action = ceil(rand*2);
            else
                [~,action] = max(Q(t-1,:));
            end
        case 'softmax'
            [softVals,action] = mySoftmax(Q(t-1,:),temp);
    end
    actions(t) = action;
    % rewards r are the per-trial EoR (ratio of reward
    % to cost)
    if (action == 1) % if chose action 1 (SR side)
        r = SR/utilityFunc2(utilityFuncValsSR(t));
        rewards(t) = SR;
        nPresses = nPresses + Ps;
        presses(t) = Ps;
        %disp(['SR trial: Ps = ' num2str(Ps) '  denom = ' num2str(utilityFunc2(utilityFuncValsSR(t)))])
    elseif (action == 2) % if chose action 2 (PR side)
        curLRtrial = curLRtrial+1;
        r = LR/utilityFunc2(utilityFuncValsPR(Pl));
        rewards(t) = LR;
        nPresses = nPresses + Pl;
        presses(t) = Pl;
        Pl = Pl + 1; % increment the PR side press cost
        %disp(['PR trial: Pl = ' num2str(Pl) '  denom = ' num2str(utilityFunc2(utilityFuncValsPR(Pl)))])
    else
        error(['action ' num2str(action) ' not recognized'])
    end
    trialEoRs(t) = r;
    %disp(['t = ' num2str(t) '  action = ' num2str(action) '  reward = ' num2str(r)])
    switch agentType
        case 'bandit'
            unchosenAction = setdiff([1 2],action);
            Q(t,action) = Q(t-1,action) + alpha*(r - Q(t-1,action));
            Q(t,unchosenAction) = Q(t-1,unchosenAction);
        case 'Qlearner'
            unchosenAction = setdiff([1 2],action);
            Q(t,action) = Q(t-1,action) + alpha*(r + gamma*max(Q(t-1,:)) - Q(t-1,action)); % update the q-table
            Q(t,unchosenAction) = Q(t-1,unchosenAction);
        otherwise
            error('agentType not recognized. How did you get here?')
    end
    %disp([num2str(Q(t,1) - Q(t-1,1)) ' ' num2str(Q(t,2) - Q(t-1,2))])
end
disp([num2str(curLRtrial) ' PR trials completed ' num2str(sum(actions==2))])
end

