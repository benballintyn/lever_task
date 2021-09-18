%analysis_120_trials
clear all;
do_save = 1;
savedir = '~/phd/lever_task/analysis_120_trials/';
if (~exist(savedir,'dir'))
    mkdir(savedir)
end
maxNtrials = 120;
mice = {'AZ04','AZ08','HP01','HP02','HP03','HP04','MA01','NS07','NS09','NS10'};
RoE_optimalities = cell(1,4);
EoR_optimalities = cell(1,4);
NL_observed = cell(1,4);
nTrials = cell(1,4);
nTrialsCompleted = cell(1,4);
nTrialsAborted = cell(1,4);
RoErandomComparison = cell(1,4);
EoRrandomComparison = cell(1,4);
totalReward = cell(1,4);
totalPresses = cell(1,4);
totalPressesSR = cell(1,4);
totalPressesLR = cell(1,4);
LRtrialsCompleted = cell(1,4);
SRtrialsCompleted = cell(1,4);
LRtrialsAborted = cell(1,4);
SRtrialsAborted = cell(1,4);
nAbortedPresses = cell(1,4);
numSRcompleted = cell(1,4);
numLRcompleted = cell(1,4);
EoRs = cell(1,4);
RoEs = cell(1,4);
instantEoR_optimalities = cell(1,4);
instantRoE_optimalities = cell(1,4);
nS = cell(1,4);
nL = cell(1,4);
nSA = cell(1,4);
nLA = cell(1,4);
for i=1:4
    nS{i} = zeros(1,maxNtrials);
    nL{i} = zeros(1,maxNtrials);
    nSA{i} = zeros(1,maxNtrials);
    nLA{i} = zeros(1,maxNtrials);
end
%RoE_NL_optimal = [70 286 448 1798]; % lever presses
RoE_NL_optimal = [70.875 286.875 448.875 1798.875]; % lever presses
%EoR_NL_optimal = [11 23 29 59]; % trials at LR
EoR_NL_optimal = [10 22 28 58];
f_NL = @(N) .5*(sqrt(8*N+9)-3);
%RoE_rstar = @(N,NL,SR,LR,Ps) (N - NL)*(SR/Ps) + LR*(-2 + sqrt(4 + 2*NL));
RoE_rstar = @(N,NL,SR,LR,Ps) (N - NL)*(SR/Ps) + LR*(.5*(sqrt(8*NL+9)-3));
%EoR_star = @(N,NL,SR,LR,Ps) (1/N)*(LR*log((NL+1)/2) + (N-NL)*(SR/Ps));
eg = double(vpa(eulergamma));
EoR_star = @(N,NL,SR,LR,Ps) (1/N)*(LR*(psi(0,NL+2) + eg - 1) + (N - NL)*(SR/Ps));

%% get random behavior optimalities
%[RoERandOptimalities,EoRRandOptimalities] = randomOptimality(1000,10000);
RoERandOptimalities = load('optimality/WT/RoERandOptimalities.mat'); RoERandOptimalities = RoERandOptimalities.RoERandOptimalities;
EoRRandOptimalities = load('optimality/WT/EoRRandOptimalities.mat'); EoRRandOptimalities = EoRRandOptimalities.EoRRandOptimalities;

%% run through processed data
count = zeros(1,4);
for i=1:length(mice)
    data = load(['processed_data/WT/' mice{i} '_ReProcessedData.mat']); data=data.ProcessedData;
    for j=1:length(data)
        if (data{j}.TotalTrialsCompleted < maxNtrials)
            continue;
        else
            sessionType = sprintf('%1$ixFR%2$i',(data{j}.Reward_L/data{j}.Reward_S),data{j}.NumPressRequired_S);
            switch sessionType
                case '2xFR6'
                    ind = 1;
                    SR = 3;
                    LR = 6;
                    Ps = 6;
                case '2xFR12'
                    ind = 2;
                    SR = 3;
                    LR = 6;
                    Ps = 12;
                case '5xFR6'
                    ind = 3;
                    SR = 3;
                    LR = 15;
                    Ps = 6;
                case '5xFR12'
                    ind = 4;
                    SR = 3;
                    LR = 15;
                    Ps = 12;
            end
            SRrewarded = find(data{j}.SideRewarded == 0); SRrewarded = SRrewarded(SRrewarded <= maxNtrials);
            LRrewarded = find(data{j}.SideRewarded == 1); LRrewarded = LRrewarded(LRrewarded <= maxNtrials);
            rewards = zeros(1,maxNtrials);
            rewards(SRrewarded) = data{j}.Reward_S;
            rewards(LRrewarded) = data{j}.Reward_L;
            EoR = zeros(1,maxNtrials);
            RoE = zeros(1,maxNtrials);
            instantEoR_optimality = zeros(1,maxNtrials);
            instantRoE_optimality = zeros(1,maxNtrials);
            for t=1:maxNtrials
                if (isnan(data{j}.SideChosen(t)))
                    continue;
                else
                    nPressesSoFar = sum(data{j}.LeverPressesByTrial(1:t));
                    EoR(t) = mean(rewards(1:t)./data{j}.LeverPressesByTrial(1:t));
                    RoE(t) = sum(rewards(1:t))/nPressesSoFar;
                    if (nPressesSoFar < RoE_NL_optimal(ind))
                        instantRoE_star = RoE_rstar(nPressesSoFar,nPressesSoFar,SR,LR,Ps)/nPressesSoFar;
                    else
                        instantRoE_star = RoE_rstar(nPressesSoFar,RoE_NL_optimal(ind),SR,LR,Ps)/nPressesSoFar;
                    end
                    if (t < EoR_NL_optimal(ind))
                        instantEoR_star = EoR_star(t,t,SR,LR,Ps);
                    else
                        instantEoR_star = EoR_star(t,EoR_NL_optimal(ind),SR,LR,Ps);
                    end
                    instantRoE_optimality(t) = RoE(t)/instantRoE_star;
                    instantEoR_optimality(t) = EoR(t)/instantEoR_star;
                    if (data{j}.SideChosen(t) == 1)
                        nL{ind}(t) = nL{ind}(t)+1;
                        if (isnan(data{j}.SideRewarded(t)))
                            nLA{ind}(t) = nLA{ind}(t)+1;
                        end
                    elseif (data{j}.SideChosen(t) == 0)
                        nS{ind}(t) = nS{ind}(t)+1;
                        if (isnan(data{j}.SideRewarded(t)))
                            nSA{ind}(t) = nSA{ind}(t)+1;
                        end
                    else
                        disp([mice{i} ' day ' num2str(j) ' trial ' num2str(t)])
                        error('side chosen not recognized')
                    end
                end
            end
            count(ind) = count(ind)+1;
            EoRs{ind}{count(ind)} = EoR;
            RoEs{ind}{count(ind)} = RoE;
            instantEoR_optimalities{ind}{count(ind)} = instantEoR_optimality;
            instantRoE_optimalities{ind}{count(ind)} = instantRoE_optimality;

            totalRobserved = sum(rewards); %data{j}.TotalRewardCollected;
            curRoE_rstar = RoE_rstar(sum(data{j}.LeverPressesByTrial(1:maxNtrials)),RoE_NL_optimal(ind),SR,LR,Ps);
            RoE_optimality = totalRobserved/curRoE_rstar;
            RoE_optimalities{ind} = [RoE_optimalities{ind} RoE_optimality];

            ok_trials = (data{j}.LeverPressesByTrial > 0);
            aborted = find(isnan(data{j}.SideRewarded)); aborted = aborted(aborted <= maxNtrials);
            presses = data{j}.LeverPressesByTrial(1:maxNtrials);
            SRchosen = find(data{j}.SideChosen == 0); SRchosen = SRchosen(SRchosen <= maxNtrials);
            LRchosen = find(data{j}.SideChosen == 1); LRchosen = LRchosen(LRchosen <= maxNtrials);
            SRcompleted = find(data{j}.SideRewarded == 0); SRcompleted = SRcompleted(SRcompleted <= maxNtrials);
            LRcompleted = find(data{j}.SideRewarded == 1); LRcompleted = LRcompleted(LRcompleted <= maxNtrials);
            SRaborted = intersect(aborted,SRchosen);
            LRaborted = intersect(aborted,LRchosen);
            SRtrialsCompleted{ind} = [SRtrialsCompleted{ind} SRcompleted];
            LRtrialsCompleted{ind} = [LRtrialsCompleted{ind} LRcompleted];
            SRtrialsAborted{ind} = [SRtrialsAborted{ind} SRaborted];
            LRtrialsAborted{ind} = [LRtrialsAborted{ind} LRaborted];
            numSRcompleted{ind} = [numSRcompleted{ind} length(SRcompleted)];
            numLRcompleted{ind} = [numLRcompleted{ind} length(LRcompleted)];
            nAbortedPresses{ind} = [nAbortedPresses{ind} sum(data{j}.LeverPressesByTrial(aborted))];
            %rewards = rewards(ok_trials);
            EoR_observed = mean(rewards./presses);
            curEoR_star = EoR_star(maxNtrials,EoR_NL_optimal(ind),SR,LR,Ps);
            EoR_optimality = EoR_observed/curEoR_star;
            if (isnan(EoR_optimality))
                disp([mice{i} ' ' num2str(j)])
            end
            EoR_optimalities{ind} = [EoR_optimalities{ind} EoR_optimality];
            NL_observed{ind} = [NL_observed{ind} length(LRcompleted)];
            nTrials{ind} = [nTrials{ind} data{j}.TotalTrialsCompleted];
            nTrialsCompleted{ind} = [nTrialsCompleted{ind} (length(LRcompleted) + length(SRcompleted))];
            nTrialsAborted{ind} = [nTrialsAborted{ind} sum(isnan(data{j}.SideRewarded(1:maxNtrials)))];
            RoErandomComparison{ind} = [RoErandomComparison{ind} RoE_optimality/mean(RoERandOptimalities{ind}(maxNtrials,:))];
            EoRrandomComparison{ind} = [EoRrandomComparison{ind} EoR_optimality/mean(EoRRandOptimalities{ind}(maxNtrials,:))];
            totalReward{ind} = [totalReward{ind} totalRobserved];
            totalPresses{ind} = [totalPresses{ind} sum(data{j}.LeverPressesByTrial(1:maxNtrials))];
            totalPressesSR{ind} = [totalPressesSR{ind} sum(data{j}.LeverPressesByTrial(SRchosen))];
            totalPressesLR{ind} = [totalPressesLR{ind} sum(data{j}.LeverPressesByTrial(LRchosen))];
        end
    end
    disp(['Done with mouse ' mice{i}])
end
for i=1:4
    pS{i} = nS{i}./(nS{i}+nL{i});
    pL{i} = nL{i}./(nS{i}+nL{i});
    pSA{i} = nSA{i}./nS{i};
    pLA{i} = nLA{i}./nL{i};
end
if (do_save)
    save([savedir 'RoE_optimalities.mat'],'RoE_optimalities','-mat')
    save([savedir 'EoR_optimalities.mat'],'EoR_optimalities','-mat')
    save([savedir 'RoEs.mat'],'RoEs','-mat')
    save([savedir 'EoRs.mat'],'EoRs','-mat')
    save([savedir 'NL_observed.mat'],'NL_observed','-mat')
    save([savedir 'nTrials.mat'],'nTrials','-mat')
    save([savedir 'nTrialsCompleted.mat'],'nTrialsCompleted','-mat')
    save([savedir 'nTrialsAborted.mat'],'nTrialsAborted','-mat')
    save([savedir 'RoERandOptimalities.mat'],'RoERandOptimalities','-mat')
    save([savedir 'EoRRandOptimalities.mat'],'EoRRandOptimalities','-mat')
    save([savedir 'RoErandomComparison.mat'],'RoErandomComparison','-mat')
    save([savedir 'EoRrandomComparison.mat'],'EoRrandomComparison','-mat')
    save([savedir 'totalReward.mat'],'totalReward','-mat')
    save([savedir 'totalPresses.mat'],'totalPresses','-mat')
    save([savedir 'totalPressesSR.mat'],'totalPressesSR','-mat')
    save([savedir 'totalPressesLR.mat'],'totalPressesLR','-mat')
    save([savedir 'pS.mat'],'pS','-mat')
    save([savedir 'pL.mat'],'pL','-mat')
    save([savedir 'pSA.mat'],'pSA','-mat')
    save([savedir 'pLA.mat'],'pLA','-mat')
end

%% Make 3D counts matrix for joint pdf objective function
nPR_counts = 15; nPR_binSize = 5;
nAborted_counts = 14; nAborted_binSize = 5;
EoR_counts = 20; EoR_binSize = .025;
mouseData3D_counts_120Trials = cell(1,4);
for i=1:4
    mouseData3D_counts_120Trials{i} = zeros(nPR_counts,nAborted_counts,EoR_counts);
end
for i=1:4
    nDataPoints = length(NL_observed{i});
    if (length(nTrialsAborted{i}) ~= nDataPoints || length(EoR_optimalities{i}) ~= nDataPoints)
        error(['Inconsistent number of data points for session type ' num2str(i)])
    end
    for j=1:nDataPoints
        % Get index for # of PR trials completed
        PR_count_ind = ceil(NL_observed{i}(j)/nPR_binSize);
        if (PR_count_ind == 0)
            PR_count_ind = 1;
        elseif (PR_count_ind > nPR_counts)
            PR_count_ind = nPR_counts;
        end
        
        % Get index for # of trials aborted
        nAborted_ind = ceil(nTrialsAborted{i}(j)/nAborted_binSize);
        if (nAborted_ind == 0)
            nAborted_ind = 1;
        elseif (nAborted_ind > nAborted_counts)
            nAborted_ind = nAborted_counts;
        end
        
        % Get index for EoR optimality
        EoR_ind = ceil(EoR_optimalities{i}(j)/EoR_binSize) - 20; % -20 is because we don't use bins for 0-.5
        if (EoR_ind == 0)
            EoR_ind = 1;
        elseif (EoR_ind > EoR_counts)
            EoR_ind = EoR_counts;
        end
        
        % add data point to joint count matrix
        mouseData3D_counts_120Trials{i}(PR_count_ind,nAborted_ind,EoR_ind) = ...
            mouseData3D_counts_120Trials{i}(PR_count_ind,nAborted_ind,EoR_ind) + 1;
    end
end
save([savedir '/mouseData3D_counts_120Trials.mat'],'mouseData3D_counts_120Trials','-mat')
    
