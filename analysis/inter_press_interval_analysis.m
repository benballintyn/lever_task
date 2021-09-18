%inter_press_interval_analysis
clear all;
animals = load('processed_data/WT/animals.mat'); animals=animals.animals;
pathToRawData = '~/phd/lever_task/raw_data/WT/';
pathToProcessedData = '~/phd/lever_task/processed_data/WT/';
sessionTypeNames = {'2xFR6','2xFR12','5xFR6','5xFR12'};
colors{1} = [0.169 0.224 0.565];
colors{2} = [0.933 0.165 0.482];
colors{3} = [0.4 0.176 0.569];
colors{4} = [0.745 0.118 0.176];
maxNpresses = 90;
for i=1:length(animals)
    ProcessedData = load([pathToProcessedData animals{i} '_ReProcessedData.mat']); ProcessedData=ProcessedData.ProcessedData;
    filenames = dir([pathToRawData animals{i} '*']);
    ipis{i} = [];
    press_positions{i} = [];
    sessionTypes{i} = [];
    isaborted{i} = [];
    isPRtrial{i} = [];
    for j=1:length(filenames)
        PR_FR_ratio = ProcessedData{j}.Reward_L / ProcessedData{j}.Reward_S;
        FR_cost = ProcessedData{j}.NumPressRequired_S;
        sessionType = [num2str(PR_FR_ratio) 'xFR' num2str(FR_cost)];
        switch sessionType
            case '2xFR6'
                sessionTypeInd = 1;
            case '2xFR12'
                sessionTypeInd = 2;
            case '5xFR6'
                sessionTypeInd = 3;
            case '5xFR12'
                sessionTypeInd = 4;
            otherwise
                error('Session Type not recognized')
        end
        rawData = load([pathToRawData filenames(j).name]); rawData=rawData.SessionData;
        disp(filenames(j).name)
        [IPI, PRtrials, Incomplete, AbortedPress] = inter_press_interval(rawData);
        for k=1:size(IPI,1)
            valid = find(~isnan(IPI(k,:)));
            ipis{i} = [ipis{i} IPI(k,valid)];
            press_positions{i} = [press_positions{i} 2:(length(valid)+1)];
            sessionTypes{i} = [sessionTypes{i} ones(1,length(valid))*sessionTypeInd];
            isaborted{i} = [isaborted{i} ones(1,length(valid))*Incomplete(k)];
            isPRtrial{i} = [isPRtrial{i} ones(1,length(valid))*PRtrials(k)];
        end
    end
end

% combine all data into vectors
allipis = [];
allpos  = [];
allSessionTypes = [];
allisaborted = [];
allisPRtrial = [];
for i=1:length(animals)
    allipis=[allipis ipis{i}];
    allpos=[allpos press_positions{i}];
    allSessionTypes = [allSessionTypes sessionTypes{i}];
    allisaborted = [allisaborted isaborted{i}];
    allisPRtrial = [allisPRtrial isPRtrial{i}];
end

% Find mean IPI, std of IPI, and number of examples for each press position
% as well as separating out presses on aborted trials
meanipis = zeros(1,maxNpresses);
nipis = zeros(1,maxNpresses);
stdipis = zeros(1,maxNpresses);
meanabortedipis = zeros(1,maxNpresses);
abortednipis = zeros(1,maxNpresses);
stdabortedipis = zeros(1,maxNpresses);
meanPRtrialipis = zeros(1,maxNpresses);
PRtrialnipis = zeros(1,maxNpresses);
stdPRtrialipis = zeros(1,maxNpresses);
meanSRtrialipis = zeros(1,maxNpresses);
SRtrialnipis = zeros(1,maxNpresses);
stdSRtrialipis = zeros(1,maxNpresses);
meanabortedPRipis = zeros(1,maxNpresses);
abortedPRnipis = zeros(1,maxNpresses);
stdabortedPRipis = zeros(1,maxNpresses);
meanabortedSRipis = zeros(1,maxNpresses);
abortedSRnipis = zeros(1,maxNpresses);
stdabortedSRipis = zeros(1,maxNpresses);
for i=2:maxNpresses
    vals = find(allpos == i);
    wasaborted = find(allisaborted(vals));
    wasPRtrial = find(allisPRtrial(vals));
    wasSRtrial = find(~allisPRtrial(vals));
    wasabortedPR = intersect(wasaborted,wasPRtrial);
    wasabortedSR = intersect(wasaborted,wasSRtrial);
    % all data
    meanipis(i) = mean(allipis(vals));
    nipis(i) = length(vals);
    stdipis(i) = std(allipis(vals));
    % aborted only trials
    meanabortedipis(i) = mean(allipis(wasaborted));
    abortednipis(i) = length(wasaborted);
    stdabortedipis(i) = std(allipis(wasaborted));
    % PR trials only
    meanPRtrialipis(i) = mean(allipis(wasPRtrial));
    PRtrialnipis(i) = length(wasPRtrial);
    stdPRtrialipis(i) = std(allipis(wasPRtrial));
    % SR trials only
    meanSRtrialipis(i) = mean(allipis(wasSRtrial));
    SRtrialnipis(i) = length(wasSRtrial);
    stdSRtrialipis(i) = std(allipis(wasSRtrial));
    % aborted PR trials only
    meanabortedPRipis(i) = mean(allipis(wasabortedPR));
    abortedPRnipis(i) = length(wasabortedPR);
    stdabortedPRipis(i) = std(allipis(wasabortedPR));
    % aborted SR trials only
    meanabortedSRipis(i) = mean(allipis(wasabortedSR));
    abortedSRnipis(i) = length(wasabortedSR);
    stdabortedSRipis(i) = std(allipis(wasabortedSR));
end
meanipis(1) = nan;
stdipis(1) = nan;
meanabortedipis(1) = nan;
stdabortedipis(1) = nan;
meanPRtrialipis(1) = nan;
stdPRtrialipis(1) = nan;
meanSRtrialipis(1) = nan;
stdSRtrialipis(1) = nan;
meanabortedPRipis(1) = nan;
stdabortedPRipis(1) = nan;
meanabortedSRipis(1) = nan;
stdabortedSRipis(1) = nan;

for i=1:4
    sessionSpecificIPIs{i} = [];
    sessionSpecificPos{i} = [];
    inds = find(allSessionTypes == i);
    sessionSpecificIPIs{i} = [sessionSpecificIPIs{i} allipis(inds)];
    sessionSpecificPos{i} = [sessionSpecificPos{i} allpos(inds)];
    sessionSpecificMeanIPIs{i} = zeros(1,maxNpresses);
    sessionSpecificStdIPIs{i} = zeros(1,maxNpresses);
    sessionSpecificNIPIs{i} = zeros(1,maxNpresses);
    for j=2:maxNpresses
        inds = find(sessionSpecificPos{i} == j);
        sessionSpecificMeanIPIs{i}(j) = mean(sessionSpecificIPIs{i}(inds));
        sessionSpecificStdIPIs{i}(j) = std(sessionSpecificIPIs{i}(inds));
        sessionSpecificNIPIs{i}(j) = length(inds);
    end
    sessionSpecificMeanIPIs{i}(1) = nan;
    sessionSpecificStdIPIs{i}(1) = nan;
    sessionSpecificNIPIs{i}(1) = nan;
end

for i=1:length(animals)
    animalSpecificMeanIPIs{i} = zeros(1,maxNpresses);
    animalSpecificStdIPIs{i} = zeros(1,maxNpresses);
    animalSpecificNIPIs{i} = zeros(1,maxNpresses);
    for j=2:maxNpresses
        inds = find(press_positions{i} == j);
        animalSpecificMeanIPIs{i}(j) = mean(ipis{i}(inds));
        animalSpecificStdIPIs{i}(j) = std(ipis{i}(inds));
        animalSpecificNIPIs{i}(j) = length(inds);
    end
    animalSpecificMeanIPIs{i}(1) = nan;
    animalSpecificStdIPIs{i}(1) = nan;
    animalSpecificNIPIs{i}(1) = nan;
end
%{
savedir = '~/phd/lever_task/IPI_analysis/';
save([savedir 'allipis.mat'],'allipis','-mat')
save([savedir 'allpos.mat'],'allpos','-mat')
save([savedir 'allSessionTypes.mat'],'allSessionTypes','-mat')
save([savedir 'allisaborted.mat'],'allisaborted','-mat')
save([savedir 'allisPRtrial.mat'],'allisPRtrial','-mat')
save([savedir 'meanipis.mat'],'meanipis','-mat')
save([savedir 'nipis.mat'],'nipis','-mat')
save([savedir 'stdipis.mat'],'stdipis','-mat')
save([savedir 'meanabortedipis.mat'],'meanabortedipis','-mat')
save([savedir 'abortednipis.mat'],'abortednipis','-mat')
save([savedir 'stdabortedipis.mat'],'stdabortedipis','-mat')
save([savedir 'meanPRtrialipis.mat'],'meanPRtrialipis','-mat')
save([savedir 'PRtrialnipis.mat'],'PRtrialnipis','-mat')
save([savedir 'stdPRtrialipis.mat'],'stdPRtrialipis','-mat')
save([savedir 'meanSRtrialipis.mat'],'meanSRtrialipis','-mat')
save([savedir 'SRtrialnipis.mat'],'SRtrialnipis','-mat')
save([savedir 'stdSRtrialipis.mat'],'stdSRtrialipis','-mat')
save([savedir 'meanabortedPRipis.mat'],'meanabortedPRipis','-mat')
save([savedir 'abortedPRnipis.mat'],'abortedPRnipis','-mat')
save([savedir 'stdabortedPRipis.mat'],'stdabortedPRipis','-mat')
save([savedir 'meanabortedSRipis.mat'],'meanabortedSRipis','-mat')
save([savedir 'abortedSRnipis.mat'],'abortedSRnipis','-mat')
save([savedir 'stdabortedSRipis.mat'],'stdabortedSRipis','-mat')
save([savedir 'sessionSpecificMeanIPIs.mat'],'sessionSpecificMeanIPIs','-mat')
save([savedir 'sessionSpecificStdIPIs.mat'],'sessionSpecificStdIPIs','-mat')
save([savedir 'sessionSpecificNIPIs.mat'],'sessionSpecificNIPIs','-mat')
save([savedir 'animalSpecificMeanIPIs.mat'],'animalSpecificMeanIPIs','-mat')
save([savedir 'animalSpecificStdIPIs.mat'],'animalSpecificStdIPIs','-mat')
save([savedir 'animalSpecificNIPIs.mat'],'animalSpecificNIPIs','-mat')

figFolder = '~/phd/lever_task/figures/WT/inter_press_interval_analysis';
figure;
subplot(1,2,1)
shadedErrorBar(1:maxNpresses,meanipis,stdipis./sqrt(nipis))
xlabel('Press position')
ylabel('Mean IPI')
title('All animals | All Sessions')
subplot(1,2,2)
plot(nipis)
xlabel('Press position')
ylabel('Total # of presses')
set(gcf,'Position',[10 10 2000 1000])
saveas(gcf,[figFolder '/AllAnimalsAllSessions.fig'],'fig')
saveas(gcf,[figFolder '/AllAnimalsAllSessions.eps'],'eps')

figure;
subplot(1,2,1);
shadedErrorBar(1:maxNpresses,meanabortedipis,stdabortedipis./sqrt(abortednipis))
xlabel('Press position')
ylabel('Mean IPI')
title('All animals | All Sessions | Aborted trials')
subplot(1,2,2)
plot(abortednipis)
xlabel('Press position')
ylabel('Total # of presses')
set(gcf,'Position',[10 10 2000 1000])
saveas(gcf,[figFolder '/AllAnimalsAllSessions_abortedTrialsOnly.fig'],'fig')
saveas(gcf,[figFolder '/AllAnimalsAllSessions_abortedTrialsOnly.eps'],'eps')

figure;
subplot(1,2,1);
shadedErrorBar(1:maxNpresses,meanPRtrialipis,stdPRtrialipis./sqrt(PRtrialnipis))
xlabel('Press position')
ylabel('Mean IPI')
title('All animals | All Sessions | PR trials')
subplot(1,2,2)
plot(PRtrialnipis)
xlabel('Press position')
ylabel('Total # of presses')
set(gcf,'Position',[10 10 2000 1000])
saveas(gcf,[figFolder '/AllAnimalsAllSessions_PRTrialsOnly.fig'],'fig')
saveas(gcf,[figFolder '/AllAnimalsAllSessions_PRTrialsOnly.eps'],'eps')

figure;
subplot(1,2,1);
shadedErrorBar(1:maxNpresses,meanSRtrialipis,stdSRtrialipis./sqrt(SRtrialnipis))
xlabel('Press position')
ylabel('Mean IPI')
title('All animals | All Sessions | SR trials')
subplot(1,2,2)
plot(SRtrialnipis)
xlabel('Press position')
ylabel('Total # of presses')
set(gcf,'Position',[10 10 2000 1000])
saveas(gcf,[figFolder '/AllAnimalsAllSessions_SRTrialsOnly.fig'],'fig')
saveas(gcf,[figFolder '/AllAnimalsAllSessions_SRTrialsOnly.eps'],'eps')

figure;
subplot(1,2,1);
shadedErrorBar(1:maxNpresses,meanabortedPRipis,stdabortedPRipis./sqrt(abortedPRnipis))
xlabel('Press position')
ylabel('Mean IPI')
title('All animals | All Sessions | aborted PR trials')
subplot(1,2,2)
plot(abortedPRnipis)
xlabel('Press position')
ylabel('Total # of presses')
set(gcf,'Position',[10 10 2000 1000])
saveas(gcf,[figFolder '/AllAnimalsAllSessions_abortedPRTrialsOnly.fig'],'fig')
saveas(gcf,[figFolder '/AllAnimalsAllSessions_abortedPRTrialsOnly.eps'],'eps')

figure;
subplot(1,2,1);
shadedErrorBar(1:maxNpresses,meanabortedSRipis,stdabortedSRipis./sqrt(abortedSRnipis))
xlabel('Press position')
ylabel('Mean IPI')
title('All animals | All Sessions | aborted SR trials')
subplot(1,2,2)
plot(abortedSRnipis)
xlabel('Press position')
ylabel('Total # of presses')
set(gcf,'Position',[10 10 2000 1000])
saveas(gcf,[figFolder '/AllAnimalsAllSessions_abortedSRTrialsOnly.fig'],'fig')
saveas(gcf,[figFolder '/AllAnimalsAllSessions_abortedSRTrialsOnly.eps'],'eps')

figure;
subplot(2,1,1)
for i=1:4
    h(i)=shadedErrorBar(1:maxNpresses,sessionSpecificMeanIPIs{i},sessionSpecificStdIPIs{i}./sqrt(sessionSpecificNIPIs{i}),'lineprops',{'markerfacecolor',colors{i}});
    hold on;
end
legend([h(1).mainLine h(2).mainLine h(3).mainLine h(4).mainLine],sessionTypeNames)
xlabel('Press position')
ylabel('Mean IPI')
title('Session specific')
subplot(2,1,2)
for i=1:4
    plot(sessionSpecificNIPIs{i},'markerfacecolor',colors{i})
    hold on;
end
legend(sessionTypeNames)
xlabel('Press position')
ylabel('Total # of presses')
set(gcf,'Position',[10 10 2000 1000])
saveas(gcf,[figFolder '/SessionSpecific.fig'],'fig')
saveas(gcf,[figFolder '/SessionSpecific.eps'],'eps')

cmap=jet;
figure;
for i=1:length(animals)
    h2(i) = shadedErrorBar(1:maxNpresses,animalSpecificMeanIPIs{i},animalSpecificStdIPIs{i}./sqrt(animalSpecificNIPIs{i}),'lineprops',{'markerfacecolor',cmap(i*floor(size(cmap,1)/length(animals)),:)});
end
legend([h2.mainLine],animals)
xlabel('Press position')
ylabel('Mean IPI')
title('Animal specific')
set(gcf,'Position',[10 10 1400 1000])
saveas(gcf,[figFolder '/AnimalSpecific.fig'],'fig')
saveas(gcf,[figFolder '/AnimalSpecific.eps'],'eps')
%}