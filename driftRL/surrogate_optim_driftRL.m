function [xmin] = surrogate_optim_driftRL(agentType,actionSelectionMethod,...
    utilityFunc1,utilityFunc2,initializationMethod,basedir,maxFunEvals)
p = inputParser;
addRequired(p,'agentType',@ischar)
addRequired(p,'actionSelectionMethod',@ischar)
addRequired(p,'utilityFunc1',@ischar)
addRequired(p,'utilityFunc2',@ischar)
addRequired(p,'initializationMethod',@ischar)
addRequired(p,'basedir',@ischar)
parse(p,agentType,actionSelectionMethod,utilityFunc1,utilityFunc2,initializationMethod,basedir)
if (~exist(basedir,'dir'))
    mkdir(basedir);
end
if (strcmp(utilityFunc1,'') && strcmp(utilityFunc2,''))
    savedir = [basedir '/' agentType '_' actionSelectionMethod '_initialization_' initializationMethod];
elseif (strcmp(utilityFunc1,'') && ~strcmp(utilityFunc2,''))
    savedir = [basedir '/' agentType '_' actionSelectionMethod '_' utilityFunc2 '_initialization_' initializationMethod];
elseif (~strcmp(utilityFunc1,'') && strcmp(utilityFunc2,''))
    savedir = [basedir '/' agentType '_' actionSelectionMethod '_' utilityFunc1 '_initialization_' initializationMethod];
else
    savedir = [basedir '/' agentType '_' actionSelectionMethod '_' utilityFunc1 '_' utilityFunc2 '_initialization_' initializationMethod];
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

runParams.agentType = agentType;
runParams.actionSelectionMethod = actionSelectionMethod;
runParams.utilityFunc1 = utilityFunc1;
runParams.utilityFunc2 = utilityFunc2;
runParams.initializationMethod = initializationMethod;
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
        curdir = [surrogateAnimalDir '/' num2str(i)];
        params = load([curdir '/params.mat']); params=params.params;
        initialPoints.X(i,:) = params;
        objective_score = load([curdir '/score.mat']); objective_score=objective_score.score;
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
switch agentType
    case 'bandit'
        lb(1) = 0; % alpha lb
        ub(1) = 1; % alpha ub
        switch actionSelectionMethod
            case 'e_greedy'
                lb(2) = 0; % epsilon lb
                ub(2) = 1; % epsilon ub
            case 'softmax'
                lb(2) = 1e-12; % temp lb
                ub(2) = 1; % temp ub
        end
        if (any(strcmp(utilityFuncs,'ansUtilityFunc')))
            lb(3) = 0; % ANS sigma lb
            ub(3) = 1; % ANS sigma ub
            lb(4) = -1; % drift rate lb
            ub(4) = 1; % drift rate ub
            lb(5) = 0; % drift noise sigma lb
            ub(5) = 1; % drift noise sigma ub
        else
            lb(3) = -1; % drift rate lb
            ub(3) = 1;  % drift rate ub
            lb(4) = 0;  % drift noise sigma lb
            ub(4) = 1;  % drift noise sigma ub
        end
        
    case 'Qlearner'
        lb(1:2) = 0; % alpha,gamma lb
        ub(1:2) = 1; % alpha,gamma ub
        switch actionSelectionMethod
            case 'e_greedy'
                lb(3) = 0; % epsilon lb
                ub(3) = 1; % epsilon ub
            case 'softmax'
                lb(3) = 1e-12; % temp lb
                ub(3) = 1; % temp ub
        end
        if (any(strcmp(utilityFuncs,'ansUtilityFunc')))
            lb(4) = 0; % ANS sigma lb
            ub(4) = 1; % ANS sigma ub
            lb(5) = -1; % drift rate lb
            ub(5) = 1; % drift rate ub
            lb(6) = 0; % drift noise sigma lb
            ub(6) = 1; % drift noise sigma ub
        else
            lb(4) = -1; % drift rate lb
            ub(4) = 1; % drift rate ub
            lb(5) = 0; % drift noise sigma lb
            ub(5) = 1; % drift noise sigma ub
        end
end

f = @(x)surrogate_optim_driftRL_inner_loop(x,savedir,agentType,actionSelectionMethod,initializationMethod,utilityFunc1,utilityFunc2);
xmin = surrogateopt(f,lb,ub,options);
end

