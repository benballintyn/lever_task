function [score] = fitFakeData_inner_loop(params,basedir,agentType,...
    actionSelectionMethod,initializationMethod,utilityFunc1,utilityFunc2,forgettingType,...
    modelType,fakeDataDir,varargin)
isValidAgentType = @(x) ismember(x,{'bandit','Qlearner'});
isValidActionSelectionMethod = @(x) ismember(x,{'e_greedy','softmax','UCB'});
isValidInitializationMethod = @(x) ismember(x,{'random','PR_biased','SR_biased','mean_reward','trained'});
isValidUtilityFunc = @(x) ismember(x,{'','ansUtilityFunc','pressUtilityFunc'});
isValidForgettingType = @(x) ismember(x,{'none','decayToInitialValues','decayToFreeParameter'});
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
addRequired(p,'modelType',isValidModelType)
addRequired(p,'fakeDataDir',@ischar)
addParameter(p,'driftType','value_based_drift',@ischar)
addParameter(p,'fullANS',false,@islogical)
addParameter(p,'noAbortANS',false,@islogical)
parse(p,params,basedir,agentType,actionSelectionMethod,initializationMethod,...
    utilityFunc1,utilityFunc2,forgettingType,modelType,fakeDataDir,varargin{:})

params = p.Results.params;

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
if (strcmp(utilityFunc1,'ansUtilityFunc'))numSR = cell(1,4); % Total # of SR trials

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

% Set up cell arrays to store behavior stats
performanceEoR = cell(1,4); % Fraction of max EoR values
nL = cell(1,4); % # of choices of PR side for each timestep
nS = cell(1,4); % # of choices of SR side for each timestep
nSaborted = cell(1,4);
nLaborted = cell(1,4);
numAborted = cell(1,4);
numLR = cell(1,4); % Total # of PR trials

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
    for nt=1:40 % Hardcoded number of sessions
        % 100 specifies the number of agents to run, 120 the number of
        % trials per session
        switch p.Results.modelType
            case 'driftRL'
                [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                    driftRLsim(sessType,120,agentType,100,...
                    actionSelectionMethod,agentParams,abortModelParams,p.Results.driftType,...
                    'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                    'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                    'Only120Trials',true,'fullANS',p.Results.fullANS);
            case 'driftRL_valueUpdate'
                [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                    driftRLsim_valueUpdate(sessType,120,agentType,100,...
                    actionSelectionMethod,agentParams,abortModelParams,p.Results.driftType,...
                    'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                    'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                    'Only120Trials',true,'fullANS',p.Results.fullANS);
            case 'logisticAbortRL'
                [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities,...
                curnSaborted,curnLaborted,curNumAborted,curPercentCompletedPR] = ...
                    logisticAbortRLsim(sessType,120,agentType,100,...
                    actionSelectionMethod,agentParams,abortModelParams,...
                    'utilityFunc1',uf1,'utilityFunc2',uf2,'initializationMethod',initializationMethod,...
                    'forgettingType',forgettingType,'forgettingParams',forgettingParams,...
                    'Only120Trials',true,'fullANS',p.Results.fullANS,'noAbortANS',p.Results.noAbortANS);
        end
        numAborted{s} = [numAborted{s} curNumAborted];
        numLR{s} = [numLR{s} curnumLR];
        performanceEoR{s} = [performanceEoR{s} EoRoptimalities];
    end
end
save([savedir '/numAborted.mat'],'numAborted','-mat')
save([savedir '/numLR.mat'],'numLR','-mat')
save([savedir '/performanceEoR.mat'],'performanceEoR','-mat')
save([savedir '/params.mat'],'params','-mat')

[score,subscores] = fake_data_objective_score(fakeDataDir,savedir);
save([savedir '/objective_score.mat'],'score','-mat')
save([savedir '/subscores.mat'],'subscores','-mat')
end

