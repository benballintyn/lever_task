function [xmin] = surrogate_policyComparison(animal,agentType,maxFunEvals)
info = load(['agent_types/' agentType '.mat']); info = info.(agentType);
info.agentParams
info.envParams
surrogateAnimalDir = [info.save_path '/' animal];
modelName = genvarname(agentType);
eval([modelName '= info'])
if (~exist(surrogateAnimalDir,'dir'))
    mkdir(surrogateAnimalDir)
    run_num = 0;
    save([surrogateAnimalDir '/run_num.mat'],'run_num')
else
    run_num = load([surrogateAnimalDir '/run_num.mat']); run_num = run_num.run_num;
end
save([surrogateAnimalDir '/' agentType '.mat'],agentType,'-mat')

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

lb = info.lower_bounds;
ub = info.upper_bounds;
f = @(x)policyComparison_inner_loop(x,info,animal);
xmin = surrogateopt(f,lb,ub,options);
end