function [score] = leverTaskPolicyComparison(animalData,animalActions,agentOpts,envOpts,agentFolder)
curAgent = leverAgent(agentOpts);
nDays = length(animalData);
actionProbs = cell(1,nDays);
allActionProbs = [];
for i=1:nDays
    curAgent.resetActionCount();
    curAgent.state = 1;
    curAgent.totalTime = 1;
    envOpts.Sreward = animalData{i}.Reward_S;
    envOpts.Lreward = animalData{i}.Reward_L;
    envOpts.nPressesForReward_S = animalData{i}.NumPressRequired_S;
    envOpts.nPressesForReward_L = animalData{i}.NumPressRequired_L(1);
    env = leverEnv(envOpts,curAgent);
    action_probabilities = zeros(1,length(animalActions{i}));
    for t=1:length(animalActions{i})
        s1 = curAgent.state;
        action = animalActions{i}(t);
        action_probabilities(t) = curAgent.getActionProb(action);
        [r,newS] = env.processAction(action);
        r = env.rewardFunc(curAgent.utilityFunc(t),r);
        curAgent.state = newS;
        curAgent.updateQ(s1,action,r,newS);
    end
    agents(i) = copy(curAgent);
    actionProbs{i} = action_probabilities;
    allActionProbs = [allActionProbs action_probabilities];
    clear action_probabilities
end
score = objective_score(allActionProbs);
save([agentFolder '/actionProbs.mat'],'actionProbs','-mat')
save([agentFolder '/agents.mat'],'agents','-mat')
save([agentFolder '/score.mat'],'score','-mat')
end

