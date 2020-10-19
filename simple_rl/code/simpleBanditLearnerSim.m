clear all;
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
epsilons = .005:.005:.1; % vector of e-greedy exploration parameters to try
%temps = .005:.005:.1; % vector of softmax exploration parameters to try
if (exist('epsilons','var')) % if using epsilon greedy method
    performance = zeros(length(alphas),length(epsilons),4); % Fraction of max RoE values
    nL = cell(length(alphas),length(epsilons),4); % # of choices of PR side for each timestep
    nS = cell(length(alphas),length(epsilons),4); % # of choices of SR side for each timestep
    pL = cell(length(alphas),length(epsilons),4); % P(PR) for each timestep. computed from nL and nS
    pS = cell(length(alphas),length(epsilons),4); % P(SR) for each timestep. computed from nL and nS
    numLR = cell(length(alphas),length(epsilons),4); % Total # of PR trials
    numSR = cell(length(alphas),length(epsilons),4); % Total # of SR trials
    nExploreParamVals = length(epsilons);
    use_e_greedy = 1;
else
    performance = zeros(length(alphas),length(temps),4);
    nL = cell(length(alphas),length(temps),4);
    nS = cell(length(alphas),length(temps),4);
    pL = cell(length(alphas),length(temps),4);
    pS = cell(length(alphas),length(temps),4);
    numLR = cell(length(alphas),length(temps),4);
    numSR = cell(length(alphas),length(temps),4);
    nExploreParamVals = length(temps);
    use_e_greedy = 0;
end

% For each parameter set, simulate nAgents q-learning agents for 1000
% trials for each session type (2xFR6,2xFR12,5xFR6,5xFR12)
for a=1:length(alphas)
    for e=1:nExploreParamVals
        parfor s=1:4
            optimalities = zeros(1,nAgents);
            % set each session types specific reward/cost parameters
            switch s
                case 1 % 2xFR6
                    SR = 3; LR = 6; Ps = 6;
                case 2 % 2xFR12
                    SR = 3; LR = 6; Ps = 12;
                case 3 % 5xFR6
                    SR = 3; LR = 15; Ps = 6;
                case 4 % 5xFR12
                    SR = 3; LR = 15; Ps = 12;
            end
            nS{a,e,s} = zeros(1,1000);
            nL{a,e,s} = zeros(1,1000);
            for i=1:nAgents
                Q = rand(1,2); % initialize a 2-element Q table with uniform random values (0,1)
                Pl = 2; % initialize PR side lever press cost
                nPresses = 0;
                rewards = zeros(1,1000);
                actions = zeros(1,1000);
                curLRtrial = 0;
                for t=1:1000
                    if (use_e_greedy) % if using epsilon greedy method
                        if (rand < epsilons(e))
                            action = ceil(rand*2); % choose action randomly
                        else
                            [~,action] = max(Q); % choose greedy action
                        end
                    else % if using the softmax method
                        [softVals,action] = mySoftmax(Q,temps(e));
                    end
                    actions(t) = action;
                    % rewards r are the per-trial EoR (ratio of reward
                    % to cost)
                    if (action == 1) % if chose action 1 (SR side)
                        r = SR/Ps;
                        rewards(t) = SR;
                        nPresses = nPresses + Ps;
                        nS{a,e,s}(t) = nS{a,e,s}(t) + 1;
                    else % if chose action 2 (PR side)
                        curLRtrial = curLRtrial+1;
                        if (curLRtrial >= 60)
                            r = LR/sum(LRPressTimeFunc2(1:(curLRtrial+1))); % denominator is a utility function fit to part of the real data from inter-press-intervals. Lever presses beyond the 60th consecutive press count for more than 1 press.
                        else
                            r = LR/Pl;
                        end
                        rewards(t) = LR;
                        nPresses = nPresses + Pl;
                        Pl = Pl + 1; % increment the PR side press cost
                        nL{a,e,s}(t) = nL{a,e,s}(t) + 1;
                    end
                    Q(action) = Q(action) + alphas(a)*(r - Q(action)); % update the q-table
                end
                curRoE_rstar = RoE_rstar(nPresses,RoE_NL_optimal(s),SR,LR,Ps);
                optimality = sum(rewards)/curRoE_rstar;
                optimalities(i) = optimality;
                numSR{a,e,s} = [numSR{a,e,s} sum(actions == 1)];
                numLR{a,e,s} = [numLR{a,e,s} sum(actions == 2)];
            end
            pS{a,e,s} = nS{a,e,s}./(nS{a,e,s} + nL{a,e,s});
            pL{a,e,s} = nL{a,e,s}./(nS{a,e,s} + nL{a,e,s});
            performance(a,e,s) = mean(optimalities);
        end
        if (use_e_greedy)
            disp(['alpha = ' num2str(alphas(a)) ' epsilon = ' num2str(epsilons(e)) ' performance = ' num2str(performance(a,e,:))])
        else
            disp(['alpha = ' num2str(alphas(a)) ' temp = ' num2str(temps(e)) ' performance = ' num2str(performance(a,e,:))])
        end
    end
end
savedir = '~/phd/lever_task/simple_rl/bandit_utility_func1_egreedy_largeSweep';
if (~exist(savedir,'dir'))
    mkdir(savedir)
end
save([savedir '/pS.mat'],'pS','-mat')
save([savedir '/pL.mat'],'pL','-mat')
save([savedir '/nS.mat'],'nS','-mat')
save([savedir '/nL.mat'],'nL','-mat')
save([savedir '/numSR.mat'],'numSR','-mat')
save([savedir '/numLR.mat'],'numLR','-mat')
save([savedir '/performance.mat'],'performance','-mat')
save([savedir '/alphas.mat'],'alphas','-mat')
if (use_e_greedy)
    save([savedir '/epsilons.mat'],'epsilons','-mat')
else
    save([savedir '/temps.mat'],'temps','-mat')
end


