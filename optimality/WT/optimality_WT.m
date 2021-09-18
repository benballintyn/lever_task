%% optimality
clear all;
do_save = 1;
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
    nS{i} = zeros(1,1000);
    nL{i} = zeros(1,1000);
    nSA{i} = zeros(1,1000);
    nLA{i} = zeros(1,1000);
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
[RoERandOptimalities,EoRRandOptimalities] = randomOptimality(1000,10000);

%% run through processed data
count = zeros(1,4);
for i=1:length(mice)
    data = load(['processed_data/' mice{i} '_ReProcessedData.mat']); data=data.ProcessedData;
    for j=1:length(data)
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
        SRrewarded = find(data{j}.SideRewarded == 0);
        LRrewarded = find(data{j}.SideRewarded == 1);
        rewards = zeros(1,data{j}.TotalTrialsCompleted);
        rewards(SRrewarded) = data{j}.Reward_S;
        rewards(LRrewarded) = data{j}.Reward_L;
        EoR = zeros(1,data{j}.TotalTrialsCompleted);
        RoE = zeros(1,data{j}.TotalTrialsCompleted);
        instantEoR_optimality = zeros(1,data{j}.TotalTrialsCompleted);
        instantRoE_optimality = zeros(1,data{j}.TotalTrialsCompleted);
        for t=1:data{j}.TotalTrialsCompleted
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
        
        totalRobserved = data{j}.TotalRewardCollected;
        curRoE_rstar = RoE_rstar(sum(data{j}.LeverPressesByTrial),RoE_NL_optimal(ind),SR,LR,Ps);
        RoE_optimality = totalRobserved/curRoE_rstar;
        RoE_optimalities{ind} = [RoE_optimalities{ind} RoE_optimality];
        
        ok_trials = (data{j}.LeverPressesByTrial > 0);
        aborted = find(isnan(data{j}.SideRewarded));
        presses = data{j}.LeverPressesByTrial(ok_trials);
        SRchosen = find(data{j}.SideChosen == 0);
        LRchosen = find(data{j}.SideChosen == 1);
        SRcompleted = find(data{j}.SideRewarded == 0);
        LRcompleted = find(data{j}.SideRewarded == 1);
        SRaborted = intersect(aborted,SRchosen);
        LRaborted = intersect(aborted,LRchosen);
        SRtrialsCompleted{ind} = [SRtrialsCompleted{ind} SRcompleted];
        LRtrialsCompleted{ind} = [LRtrialsCompleted{ind} LRcompleted];
        SRtrialsAborted{ind} = [SRtrialsAborted{ind} SRaborted];
        LRtrialsAborted{ind} = [LRtrialsAborted{ind} LRaborted];
        numSRcompleted{ind} = [numSRcompleted{ind} length(SRcompleted)];
        numLRcompleted{ind} = [numLRcompleted{ind} length(LRcompleted)];
        nAbortedPresses{ind} = [nAbortedPresses{ind} sum(data{j}.LeverPressesByTrial(aborted))];
        rewards = rewards(ok_trials);
        EoR_observed = mean(rewards./presses);
        curEoR_star = EoR_star(data{j}.TotalTrialsCompleted,EoR_NL_optimal(ind),SR,LR,Ps);
        EoR_optimality = EoR_observed/curEoR_star;
        if (isnan(EoR_optimality))
            disp([mice{i} ' ' num2str(j)])
        end
        EoR_optimalities{ind} = [EoR_optimalities{ind} EoR_optimality];
        NL_observed{ind} = [NL_observed{ind} data{j}.TotalNumTrials_LR];
        nTrials{ind} = [nTrials{ind} data{j}.TotalTrialsCompleted];
        nTrialsCompleted{ind} = [nTrialsCompleted{ind} (data{j}.TotalNumTrials_SR+data{j}.TotalNumTrials_LR)];
        nTrialsAborted{ind} = [nTrialsAborted{ind} sum(isnan(data{j}.SideRewarded))];
        RoErandomComparison{ind} = [RoErandomComparison{ind} RoE_optimality/mean(RoERandOptimalities{ind}(data{j}.TotalTrialsCompleted,:))];
        EoRrandomComparison{ind} = [EoRrandomComparison{ind} EoR_optimality/mean(EoRRandOptimalities{ind}(data{j}.TotalTrialsCompleted,:))];
        totalReward{ind} = [totalReward{ind} totalRobserved];
        totalPresses{ind} = [totalPresses{ind} sum(data{j}.LeverPressesByTrial)];
        totalPressesSR{ind} = [totalPressesSR{ind} sum(data{j}.LeverPressesByTrial(SRchosen))];
        totalPressesLR{ind} = [totalPressesLR{ind} sum(data{j}.LeverPressesByTrial(LRchosen))];
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
    save('optimality/RoE_optimalities.mat','RoE_optimalities','-mat')
    save('optimality/EoR_optimalities.mat','EoR_optimalities','-mat')
    save('optimality/RoEs.mat','RoEs','-mat')
    save('optimality/EoRs.mat','EoRs','-mat')
    save('optimality/NL_observed.mat','NL_observed','-mat')
    save('optimality/nTrials.mat','nTrials','-mat')
    save('optimality/nTrialsCompleted.mat','nTrialsCompleted','-mat')
    save('optimality/nTrialsAborted.mat','nTrialsAborted','-mat')
    save('optimality/RoERandOptimalities.mat','RoERandOptimalities','-mat')
    save('optimality/EoRRandOptimalities.mat','EoRRandOptimalities','-mat')
    save('optimality/RoErandomComparison.mat','RoErandomComparison','-mat')
    save('optimality/EoRrandomComparison.mat','EoRrandomComparison','-mat')
    save('optimality/totalReward.mat','totalReward','-mat')
    save('optimality/totalPresses.mat','totalPresses','-mat')
    save('optimality/totalPressesSR.mat','totalPressesSR','-mat')
    save('optimality/totalPressesLR.mat','totalPressesLR','-mat')
    save('optimality/pS.mat','pS','-mat')
    save('optimality/pL.mat','pL','-mat')
    save('optimality/pSA.mat','pSA','-mat')
    save('optimality/pLA.mat','pLA','-mat')
end
%% Figures
figFolder = '~/phd/lever_task/figures/optimality/';
if (~exist(figFolder,'dir'))
    mkdir(figFolder)
end
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
colors{1} = [0.169 0.224 0.565];
colors{2} = [0.933 0.165 0.482];
colors{3} = [0.4 0.176 0.569];
colors{4} = [0.745 0.118 0.176];
%%
figure;
for i=1:4
    subplot(2,2,i)
    plot(pS{i}); hold on; plot(pL{i})
    legend({'P(SR)','P(LR)'}); xlim([0 1000])
    xlabel('Trial #'); ylabel('P(Side_t = SR/LR)')
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'probLR_SR_vs_trial';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    plot(pS{i}); hold on; plot(pL{i}); plot(pSA{i}); plot(pLA{i})
    legend({'P(SR)','P(LR)','P(abort SR)','P(abort LR)'}); xlim([0 1000])
    xlabel('Trial #'); ylabel('Probability)')
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'probLR_SR_SRA_LRA_vs_trial';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    plot(pS{i}); hold on; plot(pSA{i});
    legend({'P(SR)','P(abort SR)'}); xlim([0 1000])
    xlabel('Trial #'); ylabel('Probability')
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'probSR_SRA_vs_trial';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    plot(pL{i}); hold on; plot(pLA{i});
    legend({'P(LR)','P(abort LR)'}); xlim([0 1000])
    xlabel('Trial #'); ylabel('Probability')
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'probLR_LRA_vs_trial';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    histogram(SRtrialsCompleted{i});
    xlabel('SR completed trial index'); ylabel('# of trials'); xlim([0 1000])
    title([sessionTypes{i} ' mean = ' num2str(mean(SRtrialsCompleted{i}))])
    hold on;
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'SR_trials_completed';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(SRtrialsCompleted{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('# of completed SR trials')
if (do_save)
    figname = 'SR_trials_completed_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    histogram(LRtrialsCompleted{i});
    xlabel('LR completed trial index'); ylabel('# of trials'); xlim([0 1000])
    title([sessionTypes{i} ' mean = ' num2str(mean(LRtrialsCompleted{i}))])
    hold on;
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'LR_trials_completed';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(LRtrialsCompleted{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('# of completed LR trials')
if (do_save)
    figname = 'LR_trials_completed_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    histogram(SRtrialsAborted{i});
    xlabel('SR aborted trial index'); ylabel('# of trials'); xlim([0 1000])
    title([sessionTypes{i} ' mean = ' num2str(mean(SRtrialsAborted{i}))])
    hold on;
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'SR_trials_aborted';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(SRtrialsAborted{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('SR aborted trial index')
if (do_save)
    figname = 'SR_trials_aborted_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    histogram(LRtrialsAborted{i});
    xlabel('LR aborted trial index'); ylabel('# of trials'); xlim([0 1000])
    title([sessionTypes{i} ' mean = ' num2str(mean(LRtrialsAborted{i}))])
    hold on;
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'LR_trials_aborted';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(LRtrialsAborted{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('LR aborted trial index')
if (do_save)
    figname = 'LR_trials_aborted_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

%%
figure;
for i=1:4
    subplot(2,2,i)
    for j=1:length(EoRs{i})
        plot(EoRs{i}{j})
        hold on;
    end
    ylabel('EoR(t)'); xlabel('Trial #')
    title([sessionTypes{i}])
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'EoR_vs_trial';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    for j=1:length(RoEs{i})
        plot(RoEs{i}{j})
        hold on;
    end
    ylabel('RoE(t)'); xlabel('Trial #')
    title([sessionTypes{i}])
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'RoE_vs_trial';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    for j=1:length(instantEoR_optimalities{i})
        plot(instantEoR_optimalities{i}{j})
        hold on;
    end
    ylabel('EoR(t)/EoR^*(t)')
    xlabel('Trial #')
    title(sessionTypes{i})
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'instant_EoR_optimality_by_sessionType';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    for j=1:length(instantRoE_optimalities{i})
        plot(instantRoE_optimalities{i}{j})
        hold on;
    end
    ylabel('RoE(t)/RoE^*(t)')
    xlabel('Trial #')
    title(sessionTypes{i})
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'instant_RoE_optimality_by_sessionType';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end
%%
figure;
for i=1:4
    notBoxPlot(totalPresses{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('Total # of presses')
if (do_save)
    figname = 'totalPresses_by_session_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(totalPressesSR{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('Total # of presses @SR')
if (do_save)
    figname = 'totalSRPresses_by_session_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(totalPressesLR{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('Total # of presses @LR')
if (do_save)
    figname = 'totalLRPresses_by_session_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(totalPressesLR{i}./totalPressesSR{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('LR \ SR presses ratio')
if (do_save)
    figname = 'LR_SR_press_ratio_by_session_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
maxx = max([nTrialsAborted{:}]);
for i=1:4
    subplot(2,2,i)
    linefit = polyfit(nTrialsAborted{i},RoE_optimalities{i},1);
    scatter(nTrialsAborted{i},RoE_optimalities{i},100,'.')
    hold on; x=0:1:maxx; plot(x,linefit(1)*x + linefit(2))
    xlabel('# of trials aborted'); ylabel('RoE optimality');
    text(110,.95,['y = ' num2str(linefit(1)) 'x + ' num2str(linefit(2))],'fontweight','bold');
    ylim([.5 1]);
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'nAborted_vs_RoE_by_session_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
maxx = 1000;
for i=1:4
    subplot(2,2,i)
    scatter(nTrials{i},RoE_optimalities{i},100,'.')
    hold on;
    shadedErrorBar(1:1000,mean(RoERandOptimalities{i},2),std(RoERandOptimalities{i},[],2),'lineprops','-r')
    xlabel('Total # of trials')
    ylabel('RoE optimality')
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'nTrials_vs_RoE_by_session_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(totalReward{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('Total Reward (\mul)')
if (do_save)
    figname = 'totalReward_by_sessionType_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
maxTrials = max([nTrials{:}]);
for i=1:4
    subplot(2,2,i)
    histogram(nTrials{i},'binwidth',5); xlim([0 maxTrials])
    xlabel('Total # of trials'); ylabel('# of sessions')
    title([sessionTypes{i} ' mean = ' num2str(mean(nTrials{i}))])
    set(gca,'ytick',1:max(nTrials{i}))
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'nTrials_by_sessionType';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(nTrials{i},i)
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('Total # of trials')
if (do_save)
    figname = 'nTrials_by_sessionType_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end


figure;
maxCompleted = max([nTrialsCompleted{:}]);
for i=1:4
    subplot(2,2,i)
    histogram(nTrialsCompleted{i},'binwidth',5); xlim([0 maxCompleted])
    xlabel('# of trials completed')
    title([sessionTypes{i} ' mean = ' num2str(mean(nTrialsCompleted{i}))])
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'nTrialsCompleted_by_sessionType';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(nTrialsCompleted{i},i)
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('Total # of completed trials')
if (do_save)
    figname = 'nTrialsCompleted_by_sessionType_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
maxAborted = max([nTrialsAborted{:}]);
for i=1:4
    subplot(2,2,i)
    histogram(nTrialsAborted{i},'binwidth',2); xlim([0 maxAborted])
    xlabel('# of trials aborted'); ylabel('# of sessions')
    title([sessionTypes{i} ' mean = ' num2str(mean(nTrialsAborted{i}))])
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'nTrialsAborted_by_sessionType';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(nTrialsAborted{i},i)
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('# of aborted trials')
if (do_save)
    figname = 'nTrialsAborted_by_sessionType_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end
%%
figure;
%h1=shadedErrorBar(1:1000,mean(RoERandOptimalities{1},2),std(RoERandOptimalities{1},[],2),'lineprops','-r');
h1=shadedErrorBar(1:1000,mean(RoERandOptimalities{1},2),std(RoERandOptimalities{1},[],2),'lineprops',{'color',[0.169 0.224 0.565]});
hold on;
%h2=shadedErrorBar(1:1000,mean(RoERandOptimalities{2},2),std(RoERandOptimalities{2},[],2),'lineprops','-g');
h2=shadedErrorBar(1:1000,mean(RoERandOptimalities{2},2),std(RoERandOptimalities{2},[],2),'lineprops',{'color',[0.933 0.165 0.482]});
%h3=shadedErrorBar(1:1000,mean(RoERandOptimalities{3},2),std(RoERandOptimalities{3},[],2),'lineprops','-b');
h3=shadedErrorBar(1:1000,mean(RoERandOptimalities{3},2),std(RoERandOptimalities{3},[],2),'lineprops',{'color',[0.4 0.176 0.569]});
%h4=shadedErrorBar(1:1000,mean(RoERandOptimalities{4},2),std(RoERandOptimalities{4},[],2),'lineprops','-m');
h4=shadedErrorBar(1:1000,mean(RoERandOptimalities{4},2),std(RoERandOptimalities{4},[],2),'lineprops',{'color',[0.745 0.118 0.176]});
xlabel('Total # of trials'); ylabel('Fraction of optimal reward');
legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],sessionTypes)
if (do_save)
    figname = 'randomChoiceOptimality';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end
%%
figure;
for i=1:4
    subplot(2,2,i)
    histogram(RoErandomComparison{i},'binwidth',.02);
    xlabel(['Mouse / random optimality']); ylabel('# of sessions')
    title([sessionTypes{i} ' mean = ' num2str(mean(RoErandomComparison{i}))])
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'mouse_vs_randomChoice_optimality';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(RoErandomComparison{i},i)
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('Mouse / random optimality')
if (do_save)
    figname = 'mouse_vs_randomChoice_optimality_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i);
    histogram(RoE_optimalities{i},'binwidth',.02); xlim([.5 max(1,max(RoE_optimalities{i}))])
    xlabel('Fraction of optimal reward'); ylabel('# of sessions')
    title([sessionTypes{i} ': \mu = ' num2str(mean(RoE_optimalities{i})) ' \sigma = ' num2str(std(RoE_optimalities{i}))])
end
suptitle('RoE')
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'RoE_optimality_by_sessionType';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(RoE_optimalities{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('Fraction optimal')
if (do_save)
    figname = 'RoE_optimality_by_sessionType_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i);
    histogram(EoR_optimalities{i},'binwidth',.02); xlim([.5 max(1,max(EoR_optimalities{i}))])
    xlabel('Fraction of optimal EoR'); ylabel('# of sessions')
    title([sessionTypes{i} ': \mu = ' num2str(mean(EoR_optimalities{i})) ' \sigma = ' num2str(std(EoR_optimalities{i}))])
end
suptitle('EoR')
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'EoR_optimality_by_sessionType';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(EoR_optimalities{i},i);
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('EoR optimality')
if (do_save)
    figname = 'EoR_optimality_by_sessionType_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i);
    histogram(NL_observed{i},'binwidth',2); xlabel('# of LR trials'); xlim([0 90])
    ylabel('# of sessions')
    bestNL = f_NL(RoE_NL_optimal(i));
    hold on; plot([bestNL bestNL],[0 10])
    title([sessionTypes{i} ' median = ' num2str(median(NL_observed{i})) ' optimal = ' num2str(bestNL)]);
end
set(gcf,'Position',[10 10 1500 1000])
if (do_save)
    figname = 'LR_trials_completed_by_sessionType';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end

figure;
for i=1:4
    notBoxPlot(NL_observed{i},i)
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',sessionTypes); ylabel('# of LR trials completed')
if (do_save)
    figname = 'LR_trials_completed_by_sessionType_boxPlot';
    saveas(gcf,[figFolder figname '.fig'],'fig')
    saveas(gcf,[figFolder figname '.tif'],'tiffn')
    saveas(gcf,[figFolder figname '.eps'],'eps')
end