function [xmin] = lever_task_surrogate_optim(agentType,actionSelectionMethod,...
    utilityFunc1,utilityFunc2,initializationMethod,forgettingType,scoreType,...
    modelType,basedir,maxFunEvals,varargin)
p = inputParser;
addRequired(p,'agentType',@ischar)
addRequired(p,'actionSelectionMethod',@ischar)
addRequired(p,'utilityFunc1',@ischar)
addRequired(p,'utilityFunc2',@ischar)
addRequired(p,'initializationMethod',@ischar)
addRequired(p,'forgettingType',@ischar)
addRequired(p,'scoreType',@ischar)
addRequired(p,'modelType',@ischar)
addRequired(p,'basedir',@ischar)
addRequired(p,'maxFunEvals',@isnumeric)
addParameter(p,'Only120Trials',false,@islogical)
parse(p,agentType,actionSelectionMethod,utilityFunc1,utilityFunc2,initializationMethod,...
    forgettingType,scoreType,modelType,basedir,maxFunEvals,varargin{:})

if (~exist(basedir,'dir'))
    mkdir(basedir);
end
if (strcmp(utilityFunc1,'') && strcmp(utilityFunc2,''))
    savedir = [basedir '/' agentType '_' actionSelectionMethod '_initialization_' initializationMethod '_forgettingType_' p.Results.forgettingType];
elseif (strcmp(utilityFunc1,'') && ~strcmp(utilityFunc2,''))
    savedir = [basedir '/' agentType '_' actionSelectionMethod '_' utilityFunc2 '_initialization_' initializationMethod '_forgettingType_' p.Results.forgettingType];
elseif (~strcmp(utilityFunc1,'') && strcmp(utilityFunc2,''))
    savedir = [basedir '/' agentType '_' actionSelectionMethod '_' utilityFunc1 '_initialization_' initializationMethod '_forgettingType_' p.Results.forgettingType];
else
    savedir = [basedir '/' agentType '_' actionSelectionMethod '_' utilityFunc1 '_' utilityFunc2 '_initialization_' initializationMethod '_forgettingType_' p.Results.forgettingType];
end

[status,SYSTEM_NAME] = system('hostname');
if (~status)
    switch strtrim(SYSTEM_NAME)
        case 'silmaril'
            mouseDataLoadDir = '/home/ben/phd/lever_task/cluster_code';
        case 'hpcc.brandeis.edu'
            mouseDataLoadDir = '/work/bbal/lever_task/';
    end
else
    error('System hostname not recognized')
end

runParams.agentType = p.Results.agentType;
runParams.actionSelectionMethod = p.Results.actionSelectionMethod;
runParams.utilityFunc1 = p.Results.utilityFunc1;
runParams.utilityFunc2 = p.Results.utilityFunc2;
runParams.initializationMethod = p.Results.initializationMethod;
runParams.forgettingType = p.Results.forgettingType;
runParams.scoreType = p.Results.scoreType;
runParams.savedir = savedir;
runParams.randomSeed = cputime;
rng(runParams.randomSeed)

if (~exist(runParams.savedir,'dir'))
    mkdir(runParams.savedir)
    run_num = 0;
    save([runParams.savedir '/run_num.mat'],'run_num')
else
    run_num = load([runParams.savedir '/run_num.mat']); run_num = run_num.run_num;
end

if (run_num >= 1)
    disp([num2str(run_num) ' prior runs detected'])
    for i=1:run_num
        curdir = [runParams.savedir '/' num2str(i)];
        params = load([curdir '/params.mat']); params=params.params;
        initialPoints.X(i,:) = params;
        objective_score = load([curdir '/objective_score.mat']); objective_score=objective_score.score;
        initialPoints.Fval(i) = objective_score;
    end
    options=optimoptions('surrogateopt','InitialPoints',initialPoints,'MaxFunctionEvaluations',maxFunEvals,'PlotFcn','surrogateoptplot','UseParallel',true);
    %options=optimoptions('surrogateopt','InitialPoints',initialPoints,'MaxFunctionEvaluations',maxFunEvals,'PlotFcn','surrogateoptplot','UseParallel',false);
else
    disp('No initial points')
    options=optimoptions('surrogateopt','MaxFunctionEvaluations',maxFunEvals,'PlotFcn','surrogateoptplot','UseParallel',true);
    %options=optimoptions('surrogateopt','MaxFunctionEvaluations',maxFunEvals,'PlotFcn','surrogateoptplot','UseParallel',false);
end

utilityFuncs = {utilityFunc1,utilityFunc2};
[lb,ub] = getParamBounds(p.Results.agentType,p.Results.actionSelectionMethod,...
utilityFuncs,p.Results.forgettingType,p.Results.modelType);

switch p.Results.modelType
    case 'driftRL'
    case 'driftRL_valueUpdate'
    case 'logisticAbortRL'
        f = @(x)surrogate_optim_logisticAbortRL_inner_loop(x,savedir,agentType,actionSelectionMethod,...
            initializationMethod,utilityFunc1,utilityFunc2,forgettingType,scoreType,...
            'Only120Trials',p.Results.Only120Trials);
        xmin = surrogateopt(f,lb,ub,options);
end
end

