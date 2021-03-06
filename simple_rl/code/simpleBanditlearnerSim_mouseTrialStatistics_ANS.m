clear all;
mice = {'AZ04','AZ08','HP01','HP02','HP03','HP04','MA01','NS07','NS09','NS10'};
RoE_NL_optimal = [70.875 286.875 448.875 1798.875]; % lever presses
EoR_NL_optimal = [10 22 28 58]; % trials
f_NL = @(N) .5*(sqrt(8*N+9)-3); % f(lever presses) = trials
RoE_rstar = @(N,NL,SR,LR,Ps) (N - NL)*(SR/Ps) + LR*(.5*(sqrt(8*NL+9)-3));
eg = double(vpa(eulergamma));
EoR_star = @(N,NL,SR,LR,Ps) (1/N)*(LR*(psi(0,NL+2) + eg - 1) + (N - NL)*(SR/Ps));
press_timesSR = load(['~/phd/lever_task/IPI_analysis/meanSRtrialipis.mat']); press_timesSR=press_timesSR.meanSRtrialipis;
press_timesLR = load(['~/phd/lever_task/IPI_analysis/meanPRtrialipis.mat']); press_timesLR=press_timesLR.meanPRtrialipis;
scaledPressTimesSR = press_timesSR./min(press_timesSR);
scaledPressTimesLR = press_timesLR./min(press_timesLR);
scaledPressTimesLR(end) = scaledPressTimesLR(end-1) + (scaledPressTimesLR(end-1)-scaledPressTimesLR(end-2));
scaledPressTimesLR(1) = scaledPressTimesLR(2);
LRPressTimeFit = polyfit(50:80,scaledPressTimesLR(50:80),2);
LRPressTimeFunc = @(x) LRPressTimeFit(1)*x.^2 + LRPressTimeFit(2)*x + LRPressTimeFit(3);
LRPressTimeFunc2 = @(x) (x >= 60).*LRPressTimeFunc(x) + ones(1,length(x)).*(x < 60); % Function fit to inter-press-intervals
pressUtilityFunc = @(x) sum(LRPressTimeFunc2(1:ceil(x)));
ans_sigma = .15;
ansUtilityFunc = @(x) lognrnd(log(x) - (ans_sigma^2)/2,ans_sigma);
% Create structure with necessary info on the utility functions
utilityFuncInfo.LRPressTimeFunc2 = LRPressTimeFunc2;
utilityFuncInfo.pressUtilityFunc = pressUtilityFunc;
utilityFuncInfo.ans_sigma = ans_sigma;
utilityFuncInfo.ansUtilityFunc = ansUtilityFunc;

actionSelectionMethod = 'softmax';

nAgents = 1000; % # of agents to simulate per parameter set
alphas = .05:.05:1; % vector of learning rates to try
switch actionSelectionMethod
    case 'e_greedy'
        epsilons = .005:.005:.1; % vector of e-greedy exploration parameters to try
    case 'softmax'
        temps = .01:.01:.2; % vector of softmax exploration parameters to try
    otherwise
        error('actionSelectionMethod not recognized')
end
if (exist('epsilons','var')) % if using epsilon greedy method
    performanceEoR = cell(length(alphas),length(epsilons),4); % Fraction of max EoR values
    performanceRoE = cell(length(alphas),length(epsilons),4);
    nL = cell(length(alphas),length(epsilons),4); % # of choices of PR side for each timestep
    nS = cell(length(alphas),length(epsilons),4); % # of choices of SR side for each timestep
    pL = cell(length(alphas),length(epsilons),4); % P(PR) for each timestep. computed from nL and nS
    pS = cell(length(alphas),length(epsilons),4); % P(SR) for each timestep. computed from nL and nS
    numLR = cell(length(alphas),length(epsilons),4); % Total # of PR trials
    numSR = cell(length(alphas),length(epsilons),4); % Total # of SR trials
    nExploreParamVals = length(epsilons);
    actionSelectionMethod = 'e_greedy';
else
    performanceEoR = cell(length(alphas),length(temps),4);
    performanceRoE = cell(length(alphas),length(temps),4);
    nL = cell(length(alphas),length(temps),4);
    nS = cell(length(alphas),length(temps),4);
    pL = cell(length(alphas),length(temps),4);
    pS = cell(length(alphas),length(temps),4);
    numLR = cell(length(alphas),length(temps),4);
    numSR = cell(length(alphas),length(temps),4);
    nExploreParamVals = length(temps);
    actionSelectionMethod = 'softmax';
end

mouseTrialNums = cell(4,1);
% Load mouse trial number data
for i=1:length(mice)
    data = load(['~/phd/lever_task/processed_data/WT/' mice{i} '_ReProcessedData.mat']); data=data.ProcessedData;
    for j=1:length(data)
        sessionType = sprintf('%1$ixFR%2$i',(data{j}.Reward_L/data{j}.Reward_S),data{j}.NumPressRequired_S);
        switch sessionType
            case '2xFR6'
                mouseTrialNums{1} = [mouseTrialNums{1} (data{j}.TotalNumTrials_LR+data{j}.TotalNumTrials_SR)];
            case '2xFR12'
                mouseTrialNums{2} = [mouseTrialNums{2} (data{j}.TotalNumTrials_LR+data{j}.TotalNumTrials_SR)];
            case '5xFR6'
                mouseTrialNums{3} = [mouseTrialNums{3} (data{j}.TotalNumTrials_LR+data{j}.TotalNumTrials_SR)];
            case '5xFR12'
                mouseTrialNums{4} = [mouseTrialNums{4} (data{j}.TotalNumTrials_LR+data{j}.TotalNumTrials_SR)];
        end
    end
end

% For each parameter set, simulate nAgents q-learning agents for 1000
% trials for each session type (2xFR6,2xFR12,5xFR6,5xFR12)
for a=1:length(alphas)
    for e=1:nExploreParamVals
        agentParams.alpha = alphas(a);
        switch actionSelectionMethod
            case 'e_greedy'
                agentParams.epsilon = epsilons(e);
            case 'softmax'
                agentParams.temp = temps(e);
            otherwise
                error('actionSelectionMethod not recognized')
        end
        parfor s=1:4
            % set each session types specific reward/cost parameters
            switch s
                case 1 % 2xFR6
                    sessType = '2xFR6';
                case 2 % 2xFR12
                    sessType = '2xFR12';
                case 3 % 5xFR6
                    sessType = '5xFR6';
                case 4 % 5xFR12
                    sessType = '5xFR12';
            end
            nS{a,e,s} = zeros(1,1000);
            nL{a,e,s} = zeros(1,1000);
            for nt = 1:length(mouseTrialNums{s})
                [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities] = ...
                    simpleRLsim(sessType,mouseTrialNums{s}(nt),'bandit',100,...
                    actionSelectionMethod,agentParams,'utilityFunc1',ansUtilityFunc,'utilityFunc2',pressUtilityFunc);
                nS{a,e,s} = nS{a,e,s} + curnS;
                nL{a,e,s} = nL{a,e,s} + curnL;
                numSR{a,e,s} = [numSR{a,e,s} curnumSR];
                numLR{a,e,s} = [numLR{a,e,s} curnumLR];
                performanceEoR{a,e,s} = [performanceEoR{a,e,s} EoRoptimalities];
                performanceRoE{a,e,s} = [performanceRoE{a,e,s} RoEoptimalities];
            end
            pS{a,e,s} = nS{a,e,s}./(nS{a,e,s} + nL{a,e,s});
            pL{a,e,s} = nL{a,e,s}./(nS{a,e,s} + nL{a,e,s});
            meanPerformance(s) = mean(performanceEoR{a,e,s});
        end
        switch actionSelectionMethod
            case 'e_greedy'
                disp(['alpha = ' num2str(alphas(a)) ' epsilon = ' num2str(epsilons(e)) ' performance = ' num2str(meanPerformance)])
            case 'softmax'
                disp(['alpha = ' num2str(alphas(a)) ' temp = ' num2str(temps(e)) ' performance = ' num2str(meanPerformance)])
        end
    end
end

savedir = ['~/phd/lever_task/simple_rl/bandit_' actionSelectionMethod '_largeSweep_mouseTrialStatistics_ansUtilityFunc_' num2str(ans_sigma*100) '_pressTimeFunc'];
if (~exist(savedir,'dir'))
    mkdir(savedir)
end
save([savedir '/pS.mat'],'pS','-mat')
save([savedir '/pL.mat'],'pL','-mat')
save([savedir '/nS.mat'],'nS','-mat')
save([savedir '/nL.mat'],'nL','-mat')
save([savedir '/numSR.mat'],'numSR','-mat')
save([savedir '/numLR.mat'],'numLR','-mat')
save([savedir '/performanceEoR.mat'],'performanceEoR','-mat')
save([savedir '/performanceRoE.mat'],'performanceRoE','-mat')
save([savedir '/alphas.mat'],'alphas','-mat')
save([savedir '/utilityFuncInfo.mat'],'utilityFuncInfo','-mat')
switch actionSelectionMethod
    case 'e_greedy'
        save([savedir '/epsilons.mat'],'epsilons','-mat')
    case 'softmax'
        save([savedir '/temps.mat'],'temps','-mat')
    otherwise
        error('actionSelectionMethod not recognized. Unable to save exploration parameters')
end