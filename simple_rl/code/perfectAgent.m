% perfect agent
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

%%%%%%%%%%%%%%%%%%%%%%%%% Define agent parameters %%%%%%%%%%%%%%%%%%%%%%%%%
agentParams.alpha = 1;
agentParams.gamma = 0;
agentParams.epsilon = 1;
agentParams.temp = .1;
agentParams.ans_sigma = 0;

agentType = 'Qlearner_big';
actionSelectionMethod = 'e_greedy';
valueAdjustment = 'none';
rewardType = 'divisive';

%%%%%%%%%%%%%%%%%%%%%%%%% Create utility functions %%%%%%%%%%%%%%%%%%%%%%%
% Utility function #1 (utility function from IPI analysis)
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

% Approximate number system function
ansUtilityFunc = @(x) lognrnd(log(x) - (agentParams.ans_sigma^2)/2,agentParams.ans_sigma);

%%%%%%%%%%%%%%%%%%%%%%%%% Training parameters %%%%%%%%%%%%%%%%%%%%%%%%%
nTrials = 500;
nStates = 200;
nDays   = 100;

%%%%%%%%%%%%%%%%%%%%%%%%% Do Training %%%%%%%%%%%%%%%%%%%%%%%%%
[Q,params,allSessionTypes,allActions,allRewards] = trainAgent(agentType,...
                                                           actionSelectionMethod,...
                                                           agentParams,...
                                                           'valueAdjustment',valueAdjustment,...
                                                           'rewardType',rewardType,...
                                                           'nTrials',nTrials,...
                                                           'nStates',nStates,...
                                                           'nDays',nDays);

%%%%%%%%%%%%%%%%%%%%%%%%% Test trained agent %%%%%%%%%%%%%%%%%%%%%%%%%
% Change agentParams
switch actionSelectionMethod
    case 'e_greedy'
        agentParams.epsilon = 0;
    case 'softmax'
        agentParams.temp = .000001;
end
% Set parameters for testing
nTrials = 1000;
nRepetitions = 100;
useANS = false;
valueAdjustment = 'none';
rewardType = 'divisive';
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
% Run test sessions for each sessionType
for i=1:length(sessionTypes)
    sessionType = sessionTypes{i};
    [allstates,allactions,allrewards,curNumLR] = runTrainedAgent(Q{i},agentType,actionSelectionMethod,...
                                           agentParams,nTrials,nRepetitions,sessionType,...
                                           'valueAdjustment',valueAdjustment,...
                                           'rewardType',rewardType,...
                                           'useANS',useANS);
    states{i} = allstates;
    actions{i} = allactions;
    rewards{i} = allrewards;
    numLR{i} = curNumLR;
end