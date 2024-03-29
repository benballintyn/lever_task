function [] = genFakeDataSet(SAVE_DIR,modelType,actionSelectionMethod,...
                            utilityFunc1,utilityFunc2,initializationMethod,...
                            forgettingType,varargin)
isValidActionSelectionMethod = @(x) ismember(x,{'e_greedy','softmax','UCB'});
isValidInitializationMethod = @(x) ismember(x,{'random','PR_biased','SR_biased','mean_reward','trained'});
isValidUtilityFunc = @(x) ismember(x,{'','ansUtilityFunc','pressUtilityFunc'});
isValidForgettingType = @(x) ismember(x,{'none','decayToInitialValues','decayToFreeParameter'});
isValidModelType = @(x) ismember(x,{'driftRL','driftRL_valueUpdate','logisticAbortRL'});
isValidAgentType = @(x) ismember(x,{'bandit','Qlearner'});
p = inputParser;
addRequired(p,'SAVE_DIR',@ischar)
addRequired(p,'modelType',isValidModelType)
addRequired(p,'actionSelectionMethod',isValidActionSelectionMethod)
addRequired(p,'utilityFunc1',isValidUtilityFunc)
addRequired(p,'utilityFunc2',isValidUtilityFunc)
addRequired(p,'initializationMethod',isValidInitializationMethod)
addRequired(p,'forgettingType',isValidForgettingType)
addParameter(p,'driftType','value_based_drift',@ischar)
addParameter(p,'fullANS',false,@islogical)
addParameter(p,'noAbortANS',false,@islogical)
addParameter(p,'nParamSets',10,@isnumeric)
addParameter(p,'agentType','bandit',isValidAgentType)
addParameter(p,'Only120Trials',true,@islogical)
parse(p,SAVE_DIR,modelType,actionSelectionMethod,utilityFunc1,utilityFunc2,...
    initializationMethod,forgettingType,varargin{:})

if (~exist(SAVE_DIR,'dir'))
    mkdir(SAVE_DIR)
end

utilityFuncs = {utilityFunc1,utilityFunc2};

% Loop over number of parameter sets to test
parfor i = 1:p.Results.nParamSets
    % Make directory for current parameter set runs
    CUR_SAVE_DIR = [SAVE_DIR '/' num2str(i)];
    mkdir(CUR_SAVE_DIR)
    
    % Generate parameter set
    [params,paramNames] = genParamSet(modelType,actionSelectionMethod,...
                            utilityFuncs,forgettingType);
    
    % Parse parameters into agent and abortModel param structs
    [agentParams,abortModelParams,forgettingParams] = getParamStructs(modelType,...
                                                                      params,...
                                                                      actionSelectionMethod,...
                                                                      utilityFuncs,...
                                                                      forgettingType);
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
        ansUtilityFunc = @(x,sigma) lognrnd(log(x) - (agentParams.ans_sigma^2)/2,agentParams.ans_sigma);
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
    
    %=================BEGIN SIMULATION====================================%
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
        nS{s} = zeros(1,120);
        nL{s} = zeros(1,120);
        nSaborted{s} = zeros(1,120);
        nLaborted{s} = zeros(1,120);
        for nt=1:4 % Hardcoded number of sessions
            % 100 specifies the number of agents to run, 120 the number of
            % trials per session
            switch p.Results.modelType
                case 'driftRL'
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                    curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                        driftRLsim(sessType,120,p.Results.agentType,10,...
                        actionSelectionMethod,agentParams,abortModelParams,p.Results.driftType,...
                        'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                        'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                        'Only120Trials',p.Results.Only120Trials,'fullANS',p.Results.fullANS);
                case 'driftRL_valueUpdate'
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                    curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                        driftRLsim_valueUpdate(sessType,120,p.Results.agentType,10,...
                        actionSelectionMethod,agentParams,abortModelParams,p.Results.driftType,...
                        'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                        'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                        'Only120Trials',p.Results.Only120Trials,'fullANS',p.Results.fullANS);
                case 'logisticAbortRL'
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                    curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                        logisticAbortRLsim(sessType,120,p.Results.agentType,10,...
                        actionSelectionMethod,agentParams,abortModelParams,...
                        'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                        'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                        'Only120Trials',p.Results.Only120Trials,'fullANS',p.Results.fullANS,'noAbortANS',p.Results.noAbortANS);
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
    parsave([CUR_SAVE_DIR '/nS.mat'],nS)
    parsave([CUR_SAVE_DIR '/nL.mat'],nL)
    parsave([CUR_SAVE_DIR '/nSaborted.mat'],nSaborted)
    parsave([CUR_SAVE_DIR '/nLaborted.mat'],nLaborted)
    parsave([CUR_SAVE_DIR '/numAborted.mat'],numAborted)
    parsave([CUR_SAVE_DIR '/numSR.mat'],numSR)
    parsave([CUR_SAVE_DIR '/numLR.mat'],numLR)
    parsave([CUR_SAVE_DIR '/percentCompletedPR.mat'],percentCompletedPR)
    parsave([CUR_SAVE_DIR '/performanceEoR.mat'],performanceEoR)
    parsave([CUR_SAVE_DIR '/performanceRoE.mat'],performanceRoE)
    parsave([CUR_SAVE_DIR '/pS.mat'],pS)
    parsave([CUR_SAVE_DIR '/pL.mat'],pL)
    parsave([CUR_SAVE_DIR '/params.mat'],params)
    parsave([CUR_SAVE_DIR '/paramNames.mat'],paramNames)
    
    counts3D = make3Dcounts(CUR_SAVE_DIR);
    parsave([CUR_SAVE_DIR '/counts3D.mat'],counts3D)
end

end

function [params,paramNames] = genParamSet(modelType,actionSelectionMethod,...
                            utilityFuncs,forgettingType)

[lb,ub,paramNames] = getParamBounds('bandit',actionSelectionMethod,...
utilityFuncs,forgettingType,modelType);

for i=1:length(lb)
    params(i) = rand*(ub(i) - lb(i)) + lb(i);
end
end
