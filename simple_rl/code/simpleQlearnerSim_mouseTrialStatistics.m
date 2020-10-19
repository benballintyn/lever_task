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

nAgents = 1000; % # of agents to simulate per parameter set
alphas = .05:.05:1; % vector of learning rates to try
gammas = .05:.05:1; % vector of discount factors to try
%epsilons = .005:.005:.1; % vector of e-greedy exploration parameters to try
temps = .005:.005:.1; % vector of softmax exploration parameters to try
if (exist('epsilons','var')) % if using epsilon greedy method
    performanceEoR = cell(length(alphas),length(gammas),length(epsilons),4); % Fraction of max EoR values
    performanceRoE = cell(length(alphas),length(gammas),length(epsilons),4);
    nL = cell(length(alphas),length(gammas),length(epsilons),4); % # of choices of PR side for each timestep
    nS = cell(length(alphas),length(gammas),length(epsilons),4); % # of choices of SR side for each timestep
    pL = cell(length(alphas),length(gammas),length(epsilons),4); % P(PR) for each timestep. computed from nL and nS
    pS = cell(length(alphas),length(gammas),length(epsilons),4); % P(SR) for each timestep. computed from nL and nS
    numLR = cell(length(alphas),length(gammas),length(epsilons),4); % Total # of PR trials
    numSR = cell(length(alphas),length(gammas),length(epsilons),4); % Total # of SR trials
    nExploreParamVals = length(epsilons);
    use_e_greedy = 1;
else
    performanceEoR = cell(length(alphas),length(gammas),length(temps),4);
    performanceRoE = cell(length(alphas),length(gammas),length(temps),4);
    nL = cell(length(alphas),length(gammas),length(temps),4);
    nS = cell(length(alphas),length(gammas),length(temps),4);
    pL = cell(length(alphas),length(gammas),length(temps),4);
    pS = cell(length(alphas),length(gammas),length(temps),4);
    numLR = cell(length(alphas),length(gammas),length(temps),4);
    numSR = cell(length(alphas),length(gammas),length(temps),4);
    nExploreParamVals = length(temps);
    use_e_greedy = 0;
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
    for g=1:length(gammas)
        for e=1:nExploreParamVals
            agentParams.alpha = alphas(a);
            agentParams.gamma = gammas(g);
            if (use_e_greedy)
                agentParams.epsilon = epsilons(e);
            else
                agentParams.temp = temps(e);
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
                nS{a,g,e,s} = zeros(1,1000);
                nL{a,g,e,s} = zeros(1,1000);
                for nt = 1:length(mouseTrialNums{s})
                    [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities] = ...
                        simpleRLsim(sessType,mouseTrialNums{s}(nt),'Qlearner',100,use_e_greedy,LRPressTimeFunc2,agentParams);
                    nS{a,g,e,s} = nS{a,g,e,s} + curnS;
                    nL{a,g,e,s} = nL{a,g,e,s} + curnL;
                    numSR{a,g,e,s} = [numSR{a,g,e,s} curnumSR];
                    numLR{a,g,e,s} = [numLR{a,g,e,s} curnumLR];
                    performanceEoR{a,g,e,s} = [performanceEoR{a,g,e,s} EoRoptimalities];
                    performanceRoE{a,g,e,s} = [performanceRoE{a,g,e,s} RoEoptimalities];
                end
                pS{a,g,e,s} = nS{a,g,e,s}./(nS{a,g,e,s} + nL{a,g,e,s});
                pL{a,g,e,s} = nL{a,g,e,s}./(nS{a,g,e,s} + nL{a,g,e,s});
                meanPerformance(s) = mean(performanceEoR{a,g,e,s});
            end
            if (use_e_greedy)
                disp(['alpha = ' num2str(alphas(a)) ' gamma = ' num2str(gammas(g)) ' epsilon = ' num2str(epsilons(e)) ' performance = ' num2str(meanPerformance)])
            else
                disp(['alpha = ' num2str(alphas(a)) ' gamma = ' num2str(gammas(g)) ' temp = ' num2str(temps(e)) ' performance = ' num2str(meanPerformance)])
            end
        end
    end
end
savedir = '~/phd/lever_task/simple_rl/utility_func1_softmax_largeSweep_mouseTrialStatistics';
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
save([savedir '/gammas.mat'],'gammas','-mat')
if (use_e_greedy)
    save([savedir '/epsilons.mat'],'epsilons','-mat')
else
    save([savedir '/temps.mat'],'temps','-mat')
end


