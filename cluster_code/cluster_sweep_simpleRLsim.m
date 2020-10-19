function [] = cluster_sweep_simpleRLsim(curRunID,agentType,actionSelectionMethod,utilityFunc1,utilityFunc2,datadir)
runParams.curRunID = curRunID;
runParams.agentType = agentType;
runParams.actionSelectionMethod = actionSelectionMethod;
runParams.utilityFunc1 = utilityFunc1;
runParams.utilityFunc2 = utilityFunc2;
runParams.savedir = [datadir '/' num2str(curRunID)];
runParams.randomSeed = cputime*curRunID;
rng(runParams.randomSeed)
if (~exist(runParams.savedir,'dir'))
    mkdir(runParams.savedir)
end
save([runParams.savedir '/runParams.mat'],'runParams','-mat')

% set variable ranges
alphas     = .05:.05:1;
gammas     = .05:.05:1;
epsilons   = .01:.01:.2;
temps      = .005:.005:.1;
ans_sigmas = .05:.05:.5;

% set agentParams based on curRunID and variable ranges
switch agentType
    case 'bandit'
        switch actionSelectionMethod
            case 'e_greedy'
                if (strcmp(utilityFunc1,'ansUtilityFunc') || strcmp(utilityFunc2,'ansUtilityFunc'))
                    [alpha_ind,epsilon_ind,ans_sigma_ind] = ind2sub([length(alphas) length(epsilons) length(ans_sigmas)],curRunID);
                    agentParams.alpha = alphas(alpha_ind);
                    agentParams.epsilon = epsilons(epsilon_ind);
                    agentParams.ans_sigma = ans_sigmas(ans_sigma_ind);
                else
                    [alpha_ind,epsilon_ind] = ind2sub([length(alphas) length(epsilons)],curRunID);
                    agentParams.alpha = alphas(alpha_ind);
                    agentParams.epsilon = epsilons(epsilon_ind);
                end
            case 'softmax'
                if (strcmp(utilityFunc1,'ansUtilityFunc') || strcmp(utilityFunc2,'ansUtilityFunc'))
                    [alpha_ind,temp_ind,ans_sigma_ind] = ind2sub([length(alphas) length(temps) length(ans_sigmas)],curRunID);
                    agentParams.alpha = alphas(alpha_ind);
                    agentParams.temp = temps(temp_ind);
                    agentParams.ans_sigma = ans_sigmas(ans_sigma_ind);
                else
                    [alpha_ind,temp_ind] = ind2sub([length(alphas) length(temps)],curRunID);
                    agentParams.alpha = alphas(alpha_ind);
                    agentParams.temp = temps(temp_ind);
                end
        end
    case 'Qlearner'
        switch actionSelectionMethod
            case 'e_greedy'
                if (strcmp(utilityFunc1,'ansUtilityFunc') || strcmp(utilityFunc2,'ansUtilityFunc'))
                    [alpha_ind,gamma_ind,epsilon_ind,ans_sigma_ind] = ind2sub([length(alphas) length(gammas) length(epsilons) length(ans_sigmas)],curRunID);
                    agentParams.alpha = alphas(alpha_ind);
                    agentParams.gamma = gammas(gamma_ind);
                    agentParams.epsilon = epsilons(epsilon_ind);
                    agentParams.ans_sigma = ans_sigmas(ans_sigma_ind);
                else
                    [alpha_ind,gamma_ind,epsilon_ind] = ind2sub([length(alphas) length(gammas) length(epsilons)],curRunID);
                    agentParams.alpha = alphas(alpha_ind);
                    agentParams.gamma = gammas(gamma_ind);
                    agentParams.epsilon = epsilons(epsilon_ind);
                end
            case 'softmax'
                if (strcmp(utilityFunc1,'ansUtilityFunc') || strcmp(utilityFunc2,'ansUtilityFunc'))
                    [alpha_ind,gamma_ind,temp_ind,ans_sigma_ind] = ind2sub([length(alphas) length(gammas) length(temps) length(ans_sigmas)],curRunID);
                    agentParams.alpha = alphas(alpha_ind);
                    agentParams.gamma = gammas(gamma_ind);
                    agentParams.temp = epsilons(temp_ind);
                    agentParams.ans_sigma = ans_sigmas(ans_sigma_ind);
                else
                    [alpha_ind,gamma_ind,temp_ind] = ind2sub([length(alphas) length(gammas) length(temps)],curRunID);
                    agentParams.alpha = alphas(alpha_ind);
                    agentParams.gamma = gammas(gamma_ind);
                    agentParams.temp = epsilons(temp_ind);
                end
        end
end

% Set up utility functions
%  Utility function based on press times (effort)
if (strcmp(utilityFunc1,'pressUtilityFunc') || strcmp(utilityFunc2,'pressUtilityFunc'))
    press_timesSR = load(['/home/ben/phd/lever_task/cluster_code/meanSRtrialipis.mat']); press_timesSR=press_timesSR.meanSRtrialipis;
    press_timesLR = load(['/home/ben/phd//lever_task/cluster_code/meanPRtrialipis.mat']); press_timesLR=press_timesLR.meanPRtrialipis;
    scaledPressTimesSR = press_timesSR./min(press_timesSR);
    scaledPressTimesLR = press_timesLR./min(press_timesLR);
    scaledPressTimesLR(end) = scaledPressTimesLR(end-1) + (scaledPressTimesLR(end-1)-scaledPressTimesLR(end-2));
    scaledPressTimesLR(1) = scaledPressTimesLR(2);
    LRPressTimeFit = polyfit(50:80,scaledPressTimesLR(50:80),2);
    LRPressTimeFunc = @(x) LRPressTimeFit(1)*x.^2 + LRPressTimeFit(2)*x + LRPressTimeFit(3);
    LRPressTimeFunc2 = @(x) (x >= 60).*LRPressTimeFunc(x) + ones(1,length(x)).*(x < 60); % Function fit to inter-press-intervals
    pressUtilityFunc = @(x) sum(LRPressTimeFunc2(1:ceil(x)));
end
% Approximate number system estimate function
if (strcmp(utilityFunc1,'ansUtilityFunc') || strcmp(utilityFunc2,'ansUtilityFunc'))
    ansUtilityFunc = @(x) lognrnd(log(x) - (agentParams.ans_sigma^2)/2,agentParams.ans_sigma);
end
% Set utilityFunc1 and utilityFunc2 to appropriate functions
if (strcmp(utilityFunc1,'ansUtilityFunc'))
    uf1 = ansUtilityFunc;
    agentParams.utilityFunc1 = ansUtilityFunc;
elseif (strcmp(utilityFunc1,'pressUtilityFunc'))
    uf1 = pressUtilityFunc;
    agentParams.utilityFunc1 = pressUtilityFunc;
elseif (strcmp(utilityFunc1,''))
    uf1 = @(x) x;
    agentParams.utilityFunc1 = @(x) x;
else
    error('utilityFunc1 not recognized')
end
if (strcmp(utilityFunc2,'ansUtilityFunc'))
    uf2 = ansUtilityFunc;
    agentParams.utilityFunc2 = ansUtilityFunc;
elseif (strcmp(utilityFunc2,'pressUtilityFunc'))
    uf2 = pressUtilityFunc;
    agentParams.utilityFunc2 = pressUtilityFunc;
elseif (strcmp(utilityFunc2,''))
    uf2 = @(x) x;
    agentParams.utilityFunc2 = @(x) x;
else
    error('utilityFunc2 not recognized')
end
% Save agentParams
save([runParams.savedir '/agentParams.mat'],'agentParams','-mat')

% Load mouse trial number distributions
%mouseTrialNums = load('/work/bbal/lever_task/mouseTrialNums.mat'); mouseTrialNums=mouseTrialNums.mouseTrialNums;
mouseTrialNums = load('/home/ben/phd/lever_task/cluster_code/mouseTrialNums.mat'); mouseTrialNums=mouseTrialNums.mouseTrialNums;
% Set up cell arrays to store behavior stats
performanceEoR = cell(1,4); % Fraction of max EoR values
performanceRoE = cell(1,4);
nL = cell(1,4); % # of choices of PR side for each timestep
nS = cell(1,4); % # of choices of SR side for each timestep
pL = cell(1,4); % P(PR) for each timestep. computed from nL and nS
pS = cell(1,4); % P(SR) for each timestep. computed from nL and nS
numLR = cell(1,4); % Total # of PR trials
numSR = cell(1,4); % Total # of SR trials
% For each session type and each number of trials, run simpleRLsim with
% specified parameters
for s=1:4
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
    nS{s} = zeros(1,1000);
    nL{s} = zeros(1,1000);
    for nt=1:length(mouseTrialNums{s})
        % 100 specifies the number of agents to run
        [curnumSR,curnumLR,curnS,curnL,RoEoptimalities,EoRoptimalities] = ...
                    simpleRLsim(sessType,mouseTrialNums{s}(nt),agentType,100,...
                    actionSelectionMethod,agentParams,'utilityFunc1',uf1,'utilityFunc2',uf2);
        nS{s} = nS{s} + curnS;
        nL{s} = nL{s} + curnL;
        numSR{s} = [numSR{s} curnumSR];
        numLR{s} = [numLR{s} curnumLR];
        performanceEoR{s} = [performanceEoR{s} EoRoptimalities];
        performanceRoE{s} = [performanceRoE{s} RoEoptimalities];
    end
    pS{s} = nS{s}./(nS{s} + nL{s});
    pL{s} = nL{s}./(nS{s} + nL{s});
end
save([runParams.savedir '/nS.mat'],'nS','-mat')
save([runParams.savedir '/nL.mat'],'nL','-mat')
save([runParams.savedir '/numSR.mat'],'numSR','-mat')
save([runParams.savedir '/numLR.mat'],'numLR','-mat')
save([runParams.savedir '/performanceEoR.mat'],'performanceEoR','-mat')
save([runParams.savedir '/performanceRoE.mat'],'performanceRoE','-mat')
save([runParams.savedir '/pS.mat'],'pS','-mat')
save([runParams.savedir '/pL.mat'],'pL','-mat')
end

