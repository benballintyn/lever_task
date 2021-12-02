function [numSR,numLR,nS,nL,RoEoptimalities,EoRoptimalities,nSaborted,nLaborted,...
    numAborted,percentCompletedPR] = driftRLsim(sessionType,nTrials,agentType,...
    nAgents,actionSelectionMethod,agentParams,driftParams,driftType,varargin)
RoE_NL_optimal = [70.875 286.875 448.875 1798.875]; % lever presses
EoR_NL_optimal = [10 22 28 58]; % trials
f_NL = @(N) .5*(sqrt(8*N+9)-3); % f(lever presses) = trials
RoE_rstar = @(N,NL,SR,LR,Ps) (N - NL)*(SR/Ps) + LR*(.5*(sqrt(8*NL+9)-3));
eg = double(vpa(eulergamma));
EoR_star = @(N,NL,SR,LR,Ps) (1/N)*(LR*(psi(0,NL+2) + eg - 1) + (N - NL)*(SR/Ps));

% set up inputParser object and parse inputs
positiveNoInfCheck = @(x) x > 0 && ~isinf(x);
isaFunc = @(x) isa(x,'function_handle');
p = inputParser;
addRequired(p,'sessionType',@ischar) % 2xFR6, 2xFR12, 5xFR6, 5xFR12
addRequired(p,'nTrials',positiveNoInfCheck) %  # of trials to run per agent
addRequired(p,'agentType',@ischar) % currently: 'bandit' or 'Qlearner'
addRequired(p,'nAgents',positiveNoInfCheck) % # of agents to simulate
addRequired(p,'actionSelectionMethod',@ischar) % currently: 'e_greedy' or 'softmax'
addRequired(p,'agentParams',@isstruct) % contains the parameters of the agent (alpha,gamma,epsilon/temp)
addRequired(p,'driftParams',@isstruct)
addRequired(p,'driftType',@ischar)
addParameter(p,'utilityFunc1',@(x) x,isaFunc) % optional utility function to pass lever presses through
addParameter(p,'utilityFunc2',@(x) x,isaFunc) % optional utility function to pass the results of utilityFunc1 through
addParameter(p,'initializationMethod','random',@ischar)
addParameter(p,'forgettingType','none',@ischar)
addParameter(p,'forgettingParams',[],@isnumeric)
addParameter(p,'Only120Trials',false,@islogical)
addParameter(p,'fullANS',false,@islogical)
parse(p,sessionType,nTrials,agentType,nAgents,actionSelectionMethod,agentParams,...
    driftParams,driftType,varargin{:})

if (p.Results.fullANS)
    if (~isfield(agentParams,'ans_sigma'))
        error('fullANS is selected but there is not ans_sigma parameter')
    else
        ansFunc = @(x,sigma) lognrnd(log(x) - (sigma^2)/2,sigma);
    end
end

agentType = p.Results.agentType;
actionSelectionMethod = p.Results.actionSelectionMethod;
sessionType = p.Results.sessionType;
switch p.Results.agentType
    case 'bandit'
        alpha = p.Results.agentParams.alpha;
    case 'Qlearner'
        alpha = p.Results.agentParams.alpha;
        gamma = p.Results.agentParams.gamma;
    otherwise
        error('agentType not recognized')
end
switch p.Results.actionSelectionMethod
    case 'e_greedy'
        epsilon = p.Results.agentParams.epsilon;
    case 'softmax'
        temp = p.Results.agentParams.temp;
    otherwise
        error('actionSelectionMethod not recognized')
end
switch p.Results.forgettingType
    case 'none'
        if (~isempty(p.Results.forgettingParams))
            error('forgettingType is none but forgettingParams is not empty')
        end
        useForgetting = false;
    case 'decayToInitialValues'
        if (length(p.Results.forgettingParams) ~= 1)
            error(['With forgettingType = decayToInitialValues, there should be 1 parameter but ' num2str(length(p.Results.forgettingParams)) ' were given'])
        else
            alphaF = p.Results.forgettingParams(1);
            forgetTargetCoeff = 1;
        end
        useForgetting = true;
    case 'decayToFreeParameter'
        if (length(p.Results.forgettingParams) ~= 2)
            error(['With forgettingType = decayToFreeParameter, there should be 2 parameters but ' num2str(length(p.Results.forgettingParams)) ' were given'])
        else
            alphaF = p.Results.forgettingParams(1);
            forgetTargetCoeff = p.Results.forgettingParams(2);
        end
        useForgetting = true;
    otherwise
        error('forgettingType not recognized')
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

nTrials = p.Results.nTrials;
nAgents = p.Results.nAgents;
utilityFunc1 = p.Results.utilityFunc1;
utilityFunc2 = p.Results.utilityFunc2;
initializationMethod = p.Results.initializationMethod;
numSR = [];
numLR = [];
numAborted = [];
percentCompletedPR = [];
if (p.Results.Only120Trials)
    if (nTrials ~= 120)
        error(['Only120Trials flag is true but nTrials = ' num2str(nTrials)])
    end
    nS = zeros(1,120);
    nL = zeros(1,120);
    nSaborted = zeros(1,120);
    nLaborted = zeros(1,120);
else
    nS = zeros(1,1000);
    nL = zeros(1,1000);
    nSaborted = zeros(1,1000);
    nLaborted = zeros(1,1000);
end
for i=1:nAgents
    switch initializationMethod
        case 'random'
            Q = rand(1,2); % initialize a 2-element Q table with uniform random values (0,1)
        case 'PR_biased'
            Q = [.5 .5001];
        case 'SR_biased'
            Q = [.5001 .5];
        case 'mean_reward'
            Q = [(SR + LR)/2 (SR + LR)/2];
        case 'trained'
            Q = [SR/Ps LR/2];
        otherwise
            error('initializationMethod not recognized')
    end
    initialQValues = Q;
    Pl = 2; % initialize PR side lever press cost
    nPresses = 0;
    nAborted = 0;
    PR_completed = 0;
    SR_completed = 0;
    rewards = zeros(1,nTrials);
    actions = zeros(1,nTrials);
    presses = zeros(1,nTrials);
    utilityFuncValsSR = utilityFunc1(ones(1,nTrials)*Ps);
    utilityFuncValsPR = utilityFunc1(1:(nTrials+1));
    for j=1:length(utilityFuncValsSR)
        denomsSR(j) = utilityFunc2(utilityFuncValsSR(j));
    end
    for j=1:length(utilityFuncValsPR)
        denomsPR(j) = utilityFunc2(utilityFuncValsPR(j));
    end
    for t=1:nTrials
        % Choose action based on actionSelectionMethod
        switch actionSelectionMethod
            case 'e_greedy'
                if (rand < epsilon)
                    action = ceil(rand*2);
                else
                    [~,action] = max(Q);
                end
            case 'softmax'
                [softVals,action] = mySoftmax(Q,temp);
        end
        actions(t) = action;
        % rewards r are the per-trial EoR (ratio of reward
        % to cost)
        if (action == 1) % if chose action 1 (SR side)
            switch p.Results.driftType
                case 'value_based_drift'
                    [trajectory,hitBound,abortInd] = driftProcess(driftParams.drift_rate*-(Q(1)-Q(2)),driftParams.noise_amplitude,0,10,Ps);
                case 'value_based_start'
                    [trajectory,hitBound,abortInd] = driftProcess(driftParams.drift_rate,driftParams.noise_amplitude,-(Q(1)-Q(2)),10,Ps);
            end
            if (~hitBound)
                if (p.Results.fullANS)
                    r = ansFunc(SR,agentParams.ans_sigma)/denomsSR(t);
                else
                    r = SR/denomsSR(t); %SR/utilityFunc2(utilityFuncValsSR(t));
                end
                rewards(t) = SR;
                nPresses = nPresses + Ps;
                presses(t) = Ps;
                SR_completed = SR_completed + 1;
            else
                r = 0; %abortInd*leverPressCost;
                rewards(t) = 0;
                nPresses = nPresses + abortInd;
                presses(t) = abortInd;
                nSaborted(t) = nSaborted(t) + 1;
                nAborted = nAborted + 1;
            end
            nS(t) = nS(t) + 1;
        else % if chose action 2 (PR side)
            switch p.Results.driftType
                case 'value_based_drift'
                    [trajectory,hitBound,abortInd] = driftProcess(driftParams.drift_rate*-(Q(2)-Q(1)),driftParams.noise_amplitude,0,10,Pl);
                case 'value_based_start'
                    [trajectory,hitBound,abortInd] = driftProcess(driftParams.drift_rate,driftParams.noise_amplitude,-(Q(2)-Q(1)),10,Pl);
            end
            if (~hitBound)
                if (p.Results.fullANS)
                    r = ansFunc(LR,agentParams.ans_sigma)/denomsPR(Pl);
                else
                    r = LR/denomsPR(Pl); %LR/utilityFunc2(utilityFuncValsPR(Pl));
                end
                rewards(t) = LR;
                nPresses = nPresses + Pl;
                presses(t) = Pl;
                PR_completed = PR_completed + 1;
                percentCompletedPR = [percentCompletedPR 1];
                Pl = Pl + 1; % increment the PR side press cost
            else
                r = 0; %abortInd*leverPressCost;
                rewards(t) = 0;
                nPresses = nPresses + abortInd;
                presses(t) = abortInd;
                nLaborted(t) = nLaborted(t) + 1;
                nAborted = nAborted + 1;
                percentCompletedPR = [percentCompletedPR abortInd/Pl];
            end
            nL(t) = nL(t) + 1;
        end
        unchosenAction = setdiff([1 2],action);
        switch agentType
            case 'bandit'
                Q(action) = Q(action) + alpha*(r - Q(action));
            case 'Qlearner'
                Q(action) = Q(action) + alpha*(r + gamma*max(Q) - Q(action)); % update the q-table
            otherwise
                error('agentType not recognized. How did you get here?')
        end
        if (useForgetting)
            Q(unchosenAction) = Q(unchosenAction) + alphaF*(forgetTargetCoeff*initialQValues(unchosenAction) - Q(unchosenAction));
        end
    end
    %curNumSR = sum(actions == 1);
    %curNumLR = sum(actions == 2);
    numSR = [numSR SR_completed];
    numLR = [numLR PR_completed];
    numAborted = [numAborted nAborted];
    curRoE_rstar = RoE_rstar(nPresses,RoE_NL_optimal(sessInd),SR,LR,Ps);
    curEoR_star = EoR_star(nTrials,EoR_NL_optimal(sessInd),SR,LR,Ps);
    RoEoptimality = sum(rewards)/curRoE_rstar;
    EoRoptimality = mean(rewards./presses)/curEoR_star;
    RoEoptimalities(i) = RoEoptimality;
    EoRoptimalities(i) = EoRoptimality;
end
end

