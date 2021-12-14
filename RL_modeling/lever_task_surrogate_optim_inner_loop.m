function [score] = lever_task_surrogate_optim_inner_loop(params,basedir,agentType,...
    actionSelectionMethod,initializationMethod,utilityFunc1,utilityFunc2,forgettingType,...
    scoreType,modelType,varargin)
isValidAgentType = @(x) ismember(x,{'bandit','Qlearner'});
isValidActionSelectionMethod = @(x) ismember(x,{'e_greedy','softmax'});
isValidInitializationMethod = @(x) ismember(x,{'random','PR_biased','SR_biased','mean_reward','trained'});
isValidUtilityFunc = @(x) ismember(x,{'','ansUtilityFunc','pressUtilityFunc'});
isValidForgettingType = @(x) ismember(x,{'none','decayToInitialValues','decayToFreeParameter'});
isValidScoreType = @(x) ismember(x,{'logprob_independent','logprob_joint'});
isValidModelType = @(x) ismember(x,{'driftRL','driftRL_valueUpdate','logisticAbortRL'});
p = inputParser;
addRequired(p,'params',@isnumeric)
addRequired(p,'basedir',@ischar)
addRequired(p,'agentType',isValidAgentType)
addRequired(p,'actionSelectionMethod',isValidActionSelectionMethod)
addRequired(p,'initializationMethod',isValidInitializationMethod)
addRequired(p,'utilityFunc1',isValidUtilityFunc)
addRequired(p,'utilityFunc2',isValidUtilityFunc)
addRequired(p,'forgettingType',isValidForgettingType)
addRequired(p,'scoreType',isValidScoreType)
addRequired(p,'modelType',isValidModelType)
addParameter(p,'driftType','value_based_drift',@ischar)
addParameter(p,'Only120Trials',false,@islogical)
addParameter(p,'fullANS',false,@islogical)
parse(p,params,basedir,agentType,actionSelectionMethod,initializationMethod,...
    utilityFunc1,utilityFunc2,forgettingType,scoreType,modelType,varargin{:})

params = p.Results.params;
disp(params)

[status,SYSTEM_NAME] = system('hostname');
if (~status)
    switch strtrim(SYSTEM_NAME)
        case 'silmaril'
            mouseDataLoadDir = '/home/ben/phd/lever_task/cluster_code';
        case 'hpcc.brandeis.edu'
            mouseDataLoadDir = '/work/bbal/lever_task/';
        case 'miller-lab-ubuntu2'
            mouseDataLoadDir = '/home/ben/phd/lever_task/cluster_code';
        otherwise
            error('System hostname not recognized')
    end
else
    error('system(hostname) failed')
end

runNumLoaded = 0;
while (~runNumLoaded)
    try
        run_num = load([basedir '/run_num.mat']); run_num=run_num.run_num;
        runNumLoaded = 1;
    catch e
        disp('error loading run_num')
    end
end

run_num=run_num+1; save([basedir '/run_num.mat'],'run_num')
savedir = [basedir '/' num2str(run_num)];
mkdir(savedir)

utilityFuncs = {utilityFunc1,utilityFunc2};

agentParams.alpha = params(1);
% Extract parameters from params vector
switch p.Results.agentType
    case 'bandit'
        switch actionSelectionMethod
            case 'e_greedy'
                agentParams.epsilon = params(2);
            case 'softmax'
                agentParams.temp = params(2);
        end
        if (any(strcmp(utilityFuncs,'ansUtilityFunc')))
            agentParams.ans_sigma = params(3);
            if (any(strcmp(p.Results.modelType,{'driftRL','driftRL_valueUpdate'})))
                driftParams.drift_rate = params(4);
                driftParams.noise_amplitude = params(5);
            elseif (strcmp(p.Results.modelType,'logisticAbortRL'))
                logisticParams.temp = params(4);
                logisticParams.offset = params(5);
            else
                error('modelType not recognized')
            end
            nextParamInd = 6;
        else
            if (any(strcmp(p.Results.modelType,{'driftRL','driftRL_valueUpdate'})))
                driftParams.drift_rate = params(3);
                driftParams.noise_amplitude = params(4);
            elseif (strcmp(p.Results.modelType,'logisticAbortRL'))
                logisticParams.temp = params(3);
                logisticParams.offset = params(4);
            else
                error('modelType not recognized')
            end
            nextParamInd = 5;
        end
    case 'Qlearner'
        agentParams.gamma = params(2);
        switch actionSelectionMethod
            case 'e_greedy'
                agentParams.epsilon = params(3);
            case 'softmax'
                agentParams.temp = params(3);
        end
        if (any(strcmp(utilityFuncs,'ansUtilityFunc')))
            agentParams.ans_sigma = params(4);
            if (any(strcmp(p.Results.modelType,{'driftRL','driftRL_valueUpdate'})))
                driftParams.drift_rate = params(5);
                driftParams.noise_amplitude = params(5);
            elseif (strcmp(p.Results.modelType,'logisticAbortRL'))
                logisticParams.temp = params(5);
                logisticParams.offset = params(6);
            else
                error('modelType not recognized')
            end
            nextParamInd = 7;
        else
            if (any(strcmp(p.Results.modelType,{'driftRL','driftRL_valueUpdate'})))
                driftParams.drift_rate = params(4);
                driftParams.noise_amplitude = params(5);
            elseif (strcmp(p.Results.modelType,'logisticAbortRL'))
                logisticParams.temp = params(4);
                logisticParams.offset = params(5);
            else
                error('modelType not recognized')
            end
            nextParamInd = 6;
        end
end

switch p.Results.forgettingType
    case 'none'
        if (length(params) > 6)
            error('Without forgetting, there should be at most 6 parameters')
        end
    case 'decayToInitialValues'
        forgettingParams = [params(nextParamInd)];
        if (length(params) > nextParamInd)
            error('There are extra parameters ???')
        end
    case 'decayToFreeParameter'
        forgettingParams = [params(nextParamInd:(nextParamInd+1))];
        if (length(params) > (nextParamInd+1))
            error('There are extra parameters ???')
        end
    otherwise
        error('forgettingType not recognized')
end

% Set up utility functions
%  Utility function based on press times (effort)
if (any(strcmp(utilityFuncs,'pressUtilityFunc')))
    press_timesSR = load([mouseDataLoadDir '/meanSRtrialipis.mat']); press_timesSR=press_timesSR.meanSRtrialipis;
    press_timesLR = load([mouseDataLoadDir '/meanPRtrialipis.mat']); press_timesLR=press_timesLR.meanPRtrialipis;
    scaledPressTimesSR = press_timesSR./min(press_timesSR);
    scaledPressTimesLR = press_timesLR./min(press_timesLR);
    scaledPressTimesLR(end) = scaledPressTimesLR(end-1) + (scaledPressTimesLR(end-1)-scaledPressTimesLR(end-2));
    scaledPressTimesLR(1) = scaledPressTimesLR(2);
    LRPressTimeFit = polyfit(50:80,scaledPressTimesLR(50:80),2);
    LRPressTimeFunc = @(x) LRPressTimeFit(1)*x.^2 + LRPressTimeFit(2)*x + LRPressTimeFit(3);
    LRPressTimeFunc2 = @(x) (x >= 60).*LRPressTimeFunc(x) + ones(1,length(x)).*(x < 60); % Function fit to inter-press-intervals
    pressUtilityFunc = @(x) sum(LRPressTimeFunc2(1:ceil(x)));
end
% Approximate number system estimate function
if (any(strcmp(utilityFuncs,'ansUtilityFunc')))
    ansUtilityFunc = @(x,sigma) lognrnd(log(x) - (agentParams.sigma^2)/2,agentParams.sigma);
end
% Set utilityFunc1 and utilityFunc2 to appropriate functions
if (strcmp(utilityFunc1,'ansUtilityFunc'))
    uf1 = ansUtilityFunc;
elseif (strcmp(utilityFunc1,'pressUtilityFunc'))
    uf1 = pressUtilityFunc;
elseif (strcmp(utilityFunc1,''))
    uf1 = @(x) x;
else
    error('utilityFunc1 not recognized')
end
if (strcmp(utilityFunc2,'ansUtilityFunc'))
    uf2 = ansUtilityFunc;
elseif (strcmp(utilityFunc2,'pressUtilityFunc'))
    uf2 = pressUtilityFunc;
elseif (strcmp(utilityFunc2,''))
    uf2 = @(x) x;
else
    error('utilityFunc2 not recognized')
end

if (~p.Results.Only120Trials)
    % Load mouse trial number distributions
    mouseTrialNums = load([mouseDataLoadDir '/mouseTrialNums.mat']); mouseTrialNums=mouseTrialNums.mouseTrialNums;
end

% Set up cell arrays to store behavior stats
performanceEoR = cell(1,4); % Fraction of max EoR values
performanceRoE = cell(1,4);
nL = cell(1,4); % # of choices of PR side for each timestep
nS = cell(1,4); % # of choices of SR side for each timestep
nSaborted = cell(1,4);
nLaborted = cell(1,4);
numAborted = cell(1,4);
percentCompletedPR = cell(1,4);
pL = cell(1,4); % P(PR) for each timestep. computed from nL and nS
pS = cell(1,4); % P(SR) for each timestep. computed from nL and nS
numLR = cell(1,4); % Total # of PR trials
numSR = cell(1,4); % Total # of SR trials

% For each session type and each number of trials, run corresponding sim code with
% specified parameters
for s=1:4
    switch s
        case 1 % 2xFR6
            sessType = '2xFR6';
        case 2 % 2xFR12
            sessType = '2xFR12';
        case 3 % 5xFR6
            sessType = '5xFR6';
        case 4 % 5xFR12
            sessType = '5xFR12';
    end
    if (p.Results.Only120Trials) % Simulate only 120 trials
        nS{s} = zeros(1,120);
        nL{s} = zeros(1,120);
        nSaborted{s} = zeros(1,120);
        nLaborted{s} = zeros(1,120);
        for nt=1:40 % Hardcoded number of sessions
            % 100 specifies the number of agents to run, 120 the number of
            % trials per session
            switch p.Results.modelType
                case 'driftRL'
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                    curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                        driftRLsim(sessType,120,agentType,100,...
                        actionSelectionMethod,agentParams,driftParams,p.Results.driftType,...
                        'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                        'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                        'Only120Trials',p.Results.Only120Trials,'fullANS',p.Results.fullANS);
                case 'driftRL_valueUpdate'
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                    curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                        driftRLsim_valueUpdate(sessType,120,agentType,100,...
                        actionSelectionMethod,agentParams,driftParams,p.Results.driftType,...
                        'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                        'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                        'Only120Trials',p.Results.Only120Trials,'fullANS',p.Results.fullANS);
                case 'logisticAbortRL'
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                    curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                        logisticAbortRLsim(sessType,120,agentType,100,...
                        actionSelectionMethod,agentParams,logisticParams,...
                        'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                        'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                        'Only120Trials',p.Results.Only120Trials,'fullANS',p.Results.fullANS);
            end
            nS{s} = nS{s} + curnS;
            nL{s} = nL{s} + curnL;
            nSaborted{s} = nSaborted{s} + curnSaborted;
            nLaborted{s} = nLaborted{s} + curnLaborted;
            numAborted{s} = [numAborted{s} curNumAborted];
            numSR{s} = [numSR{s} curnumSR];
            numLR{s} = [numLR{s} curnumLR];
            percentCompletedPR{s} = [percentCompletedPR{s} curPercentCompletedPR];
            performanceEoR{s} = [performanceEoR{s} EoRoptimalities];
            performanceRoE{s} = [performanceRoE{s} RoEoptimalities];
        end
        pS{s} = nS{s}./(nS{s} + nL{s});
        pL{s} = nL{s}./(nS{s} + nL{s});
    else % Simulate same number of trials as total that mice completed
        nS{s} = zeros(1,1000);
        nL{s} = zeros(1,1000);
        nSaborted{s} = zeros(1,1000);
        nLaborted{s} = zeros(1,1000);
        for nt=1:length(mouseTrialNums{s})
            % 100 specifies the number of agents to run
            switch p.Results.modelType
                case 'driftRL'
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                    curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                        driftRLsim(sessType,mouseTrialNums{s}(nt),agentType,100,...
                        actionSelectionMethod,agentParams,driftParams,p.Results.driftType,...
                        'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                        'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                        'Only120Trials',p.Results.Only120Trials,'fullANS',p.Results.fullANS);
                case 'driftRL_valueUpdate'
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                    curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                        driftRLsim_valueUpdate(sessType,mouseTrialNums{s}(nt),agentType,100,...
                        actionSelectionMethod,agentParams,driftParams,p.Results.driftType,...
                        'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                        'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                        'Only120Trials',p.Results.Only120Trials,'fullANS',p.Results.fullANS);
                case 'logisticAbortRL'
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                    curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                        logisticAbortRLsim(sessType,mouseTrialNums{s}(nt),agentType,100,...
                        actionSelectionMethod,agentParams,logisticParams,...
                        'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                        'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                        'Only120Trials',p.Results.Only120Trials,'fullANS',p.Results.fullANS);
            end
            nS{s} = nS{s} + curnS;
            nL{s} = nL{s} + curnL;
            nSaborted{s} = nSaborted{s} + curnSaborted;
            nLaborted{s} = nLaborted{s} + curnLaborted;
            numAborted{s} = [numAborted{s} curNumAborted];
            numSR{s} = [numSR{s} curnumSR];
            numLR{s} = [numLR{s} curnumLR];
            percentCompletedPR{s} = [percentCompletedPR{s} curPercentCompletedPR];
            performanceEoR{s} = [performanceEoR{s} EoRoptimalities];
            performanceRoE{s} = [performanceRoE{s} RoEoptimalities];
        end
        pS{s} = nS{s}./(nS{s} + nL{s});
        pL{s} = nL{s}./(nS{s} + nL{s});
    end
    
end
save([savedir '/nS.mat'],'nS','-mat')
save([savedir '/nL.mat'],'nL','-mat')
save([savedir '/nSaborted.mat'],'nSaborted','-mat')
save([savedir '/nLaborted.mat'],'nLaborted','-mat')
save([savedir '/numAborted.mat'],'numAborted','-mat')
save([savedir '/numSR.mat'],'numSR','-mat')
save([savedir '/numLR.mat'],'numLR','-mat')
save([savedir '/percentCompletedPR.mat'],'percentCompletedPR','-mat')
save([savedir '/performanceEoR.mat'],'performanceEoR','-mat')
save([savedir '/performanceRoE.mat'],'performanceRoE','-mat')
save([savedir '/pS.mat'],'pS','-mat')
save([savedir '/pL.mat'],'pL','-mat')
save([savedir '/params.mat'],'params','-mat')

[score,subscores] = rl_objective_score(savedir,scoreType,'Only120Trials',p.Results.Only120Trials);
save([savedir '/objective_score.mat'],'score','-mat')
save([savedir '/subscores.mat'],'subscores','-mat')
end

