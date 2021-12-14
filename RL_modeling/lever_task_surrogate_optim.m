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
addParameter(p,'driftType','value_based_drift',@ischar)
addParameter(p,'Only120Trials',false,@islogical)
addParameter(p,'fullANS',false,@islogical)
parse(p,agentType,actionSelectionMethod,utilityFunc1,utilityFunc2,initializationMethod,...
    forgettingType,scoreType,modelType,basedir,maxFunEvals,varargin{:})

if (~any(strcmp({p.Results.utilityFunc1,p.Results.utilityFunc2},'ansUtilityFunc')) && p.Results.fullANS)
    error('If fullANS == true, one of the utility functions must be ansUtilityFunc')
end

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
if (p.Results.fullANS)
    savedir = [savedir '_fullANS'];
end
savedir = [savedir '/optimized/'];

runParams.agentType = p.Results.agentType;
runParams.actionSelectionMethod = p.Results.actionSelectionMethod;
runParams.utilityFunc1 = p.Results.utilityFunc1;
runParams.utilityFunc2 = p.Results.utilityFunc2;
runParams.initializationMethod = p.Results.initializationMethod;
runParams.forgettingType = p.Results.forgettingType;
runParams.scoreType = p.Results.scoreType;
runParams.savedir = savedir;
runParams.randomSeed = cputime;
if (p.Results.fullANS)
    runParams.fullANS = true;
else
    runParams.fullANS = false;
end
rng(runParams.randomSeed)

utilityFuncs = {utilityFunc1,utilityFunc2};
[lb,ub,paramNames] = getParamBounds(p.Results.agentType,p.Results.actionSelectionMethod,...
utilityFuncs,p.Results.forgettingType,p.Results.modelType);

runParams.paramNames = paramNames;
runParams.lower_bounds = lb;
runParams.upper_bounds = ub;

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

f = @(x) lever_task_surrogate_optim_inner_loop(x,savedir,agentType,actionSelectionMethod,...
            initializationMethod,utilityFunc1,utilityFunc2,forgettingType,scoreType,p.Results.modelType,...
            'Only120Trials',p.Results.Only120Trials,'driftType',p.Results.driftType,'fullANS',p.Results.fullANS);
xmin = surrogateopt(f,lb,ub,options);
end

