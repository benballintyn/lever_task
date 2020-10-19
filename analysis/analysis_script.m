% analysis_script
animals = load('processed_data/animals.mat'); animals=animals.animals;
driveDir = '/media/ben/Manwe/phd/lever_task/agent_data';
agentTypes = {'e_greedy_noModel','e_greedy_UCB_noModel','softmax_noModel','softmax_UCB_noModel'};
agentTypes2 = {'\epsilon-greedy','\epsilon-greedy UCB','softmax','softmax UCB'};

topScores=cell(1,length(agentTypes));
for i=1:length(agentTypes)
    curAgentInfo = load(['agent_types/' agentTypes{i} '.mat']); curAgentInfo=curAgentInfo.(agentTypes{i});
    disp(agentTypes{i})
    bestVals = [];
    clear vals;
    for j=1:length(animals)
        [paramVals,paramNames,scores,subscores] = getParamVals(driveDir,agentTypes{i},animals{j});
        [~,inds] = sort(scores);
        topScores{i} = [topScores{i} scores(inds(1:10))];
        bestVals = [bestVals; paramVals(inds(1:10),:)];
    end
    figure;
    for j=1:size(bestVals,2)
        subplot(1,size(bestVals,2),j)
        scatter(ones(1,size(bestVals,1))*j,bestVals(:,j),'.')
        set(gca,'xtick',j,'xticklabels',paramNames{j},'xticklabelrotation',20)
        ylim([curAgentInfo.lower_bounds(j) curAgentInfo.upper_bounds(j)])
    end
    suptitle(agentTypes2{i})
end
figure;
violin(topScores)
set(gca,'xtick',1:length(agentTypes),'xticklabels',agentTypes,'xticklabelrotation',20)