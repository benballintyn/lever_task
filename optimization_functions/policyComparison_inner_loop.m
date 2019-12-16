function [score] = policyComparison_inner_loop(params,info,animal)
paramNames = [info.agentParams info.envParams];
params
savedir = [info.save_path '/' animal];
runNumLoaded = 0;
while (~runNumLoaded)
    try
        run_num = load([savedir '/run_num.mat']); run_num=run_num.run_num;
        runNumLoaded = 1;
    catch e
        disp('error loading run_num')
    end
end
run_num=run_num+1; save([savedir '/run_num.mat'],'run_num')
rundir = [savedir '/' num2str(run_num)]; mkdir(rundir)
for i=1:length(info.agentParams)
    agentOpts.(paramNames{i}) = params(i);
end
agentOpts.actionSelectionMethod = info.actionSelectionMethod;
agentOpts.updateMethod = info.updateMethod;
%agentOpts.nUtilityStates = info.nUtilityStates;
%agentOpts.nFatigueStates = info.nFatigueStates;
if (isfield(info,'valueAdjustment'))
    agentOpts.valueAdjustment = info.valueAdjustment;
end

for i=1:length(info.envParams)
    curInd = length(info.agentParams)+i;
    envOpts.(paramNames{curInd}) = params(curInd);
end
envOpts.rewardFunc = info.rewardFunc;
if (isfield(info,'itiCost'))
    envOpts.itiCost = info.itiCost;
end
if (isfield(info,'leverPressCost'))
    envOpts.leverPressCost = info.leverPressCost;
end
if (isfield(info,'rewardFunc'))
    envOpts.rewardFunc = info.rewardFunc;
else
    error('No reward function specified')
end
%envOpts.increaseUtilityProb = info.increaseUtilityProb;
%envOpts.decreaseUtilityProb = info.decreaseUtilityProb;
%envOpts.increaseFatigueProb = info.increaseFatigueProb;
%envOpts.decreaseFatigueProb = info.decreaseFatigueProb;
animalData = load(['processed_data/' animal '_ReProcessedData.mat']); animalData=animalData.ProcessedData;
animalActions = load(['processed_data/' animal '_ActionSequences.mat']); animalActions=animalActions.ActionSequences;
score = leverTaskPolicyComparison(animalData,animalActions,agentOpts,envOpts,rundir);

disp(['score = ' num2str(score)])
save([rundir '/params.mat'],'params','-mat')
save([rundir '/agentOpts.mat'],'agentOpts','-mat')
save([rundir '/envOpts.mat'],'envOpts','-mat')
end