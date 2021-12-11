function [scores] = grid_search(agentType,actionSelectionMethod,...
    utilityFunc1,utilityFunc2,initializationMethod,forgettingType,scoreType,...
    modelType,basedir,varargin)
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
addParameter(p,'driftType','value_based_drift',@ischar)
addParameter(p,'Only120Trials',false,@islogical)
addParameter(p,'fullANS',false,@islogical)
parse(p,agentType,actionSelectionMethod,utilityFunc1,utilityFunc2,initializationMethod,...
    forgettingType,scoreType,modelType,basedir,varargin{:})

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
savedir = [savedir '/grid_search/'];

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
runParams.modelType = p.Results.modelType;
runParams.Only120Trials = p.Results.Only120Trials;
runParams.driftType = p.Results.driftType;
runParams.savedir = savedir;
runParams.randomSeed = cputime;
rng(runParams.randomSeed)

if (~exist(runParams.savedir,'dir'))
    mkdir(runParams.savedir)
end
save([runParams.savedir '/runParams.mat'],'runParams','-mat')

utilityFuncs = {utilityFunc1,utilityFunc2};
[lb,ub] = getParamBounds(p.Results.agentType,p.Results.actionSelectionMethod,...
utilityFuncs,p.Results.forgettingType,p.Results.modelType);

N_PARAMS = length(lb);
N_SAMPLES_PER_PARAM = 4;

siz = ones(1,N_PARAMS)*N_SAMPLES_PER_PARAM;


agentType = p.Results.agentType;
actionSelectionMethod = p.Results.actionSelectionMethod;
initializationMethod = p.Results.initializationMethod;
utilityFunc1 = p.Results.utilityFunc1;
utilityFunc2 = p.Results.utilityFunc2;
forgettingType = p.Results.forgettingType;
scoreType = p.Results.scoreType;
modelType = p.Results.modelType;
Only120Trials = p.Results.Only120Trials;
driftType = p.Results.driftType;

parfor i=1:N_SAMPLES_PER_PARAM^N_PARAMS
    if (exist([savedir '/' num2str(i)],'dir'))
        continue;
    end
    tempInds = cell(size(siz));
    [tempInds{:}] = ind2sub(siz,i);
    paramInds = [tempInds{:}];
    params = getParams(lb,ub,paramInds,N_PARAMS,N_SAMPLES_PER_PARAM);
    
    if (any(isnan(params)))
        error('A parameter was assigned NaN')
    end
    
    [score] = run_param_set(params,savedir,agentType,actionSelectionMethod,...
        initializationMethod,utilityFunc1,utilityFunc2,forgettingType,scoreType,...
        modelType,i,'Only120Trials',Only120Trials,'driftType',driftType,'fullANS',p.Results.fullANS);
    
    scores(i) = score;
end
end

function params = getParams(lb,ub,paramInds,nparams,nsamplesperparam)
    params = nan(1,nparams);
    for i=1:nparams
        params(i) = ((paramInds(i) - 1)/(nsamplesperparam - 1))*(ub(i) - lb(i)) + lb(i);
    end
end
