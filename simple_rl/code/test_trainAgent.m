% train_agents
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

% Set agentParams
agentParams.alpha = .05;
agentParams.gamma = 0;
agentParams.epsilon = 1;
agentParams.ans_sigma = .5;

% set up utility functions
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

ansUtilityFunc = @(x) lognrnd(log(x) - (agentParams.ans_sigma^2)/2,agentParams.ans_sigma);

% Train agent
[Q,params,sessionTypes,allActions,allRewards] = trainAgent('Qlearner_big',...
                                                           'e_greedy',...
                                                           agentParams,...
                                                           'utilityFunc1',ansUtilityFunc,...
                                                           'utilityFunc2',pressUtilityFunc,...
                                                           'rewardType','divisive',...
                                                           'nDays',120);