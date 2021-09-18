%%
clear all; close all
saveDir = '~/phd/lever_task/driftRL/abort_analysis/';
figDir = [saveDir 'figures/'];
%%
mice = load('processed_data/WT/animals.mat'); mice=mice.animals;
abort_countPR_all = zeros(100,100);
total_countPR_all = zeros(100,100);
abort_countFR_all = zeros(2,12);
total_countFR_all = zeros(2,12);

abort_countPR = cell(1,4);
total_countPR = cell(1,4);
abort_countFR = cell(1,4);
total_countFR = cell(1,4);

percentCompletedPR = cell(1,4);
percentCompletedFR = cell(1,4);
nAborted = cell(1,4);
dataMatrix = [];
for i=1:4
    abort_countPR{i} = zeros(100,100);
    total_countPR{i} = zeros(100,100);
    if (i == 1 || i == 3)
        abort_countFR{i} = zeros(1,6);
        total_countFR{i} = zeros(1,6);
    elseif (i == 2 || i == 4)
        abort_countFR{i} = zeros(1,12);
        total_countFR{i} = zeros(1,12);
    end
end

count = 0;
mouseNames = {};
for i=1:length(mice)
    data = load(['processed_data/WT/' mice{i} '_ReProcessedData.mat']); data=data.ProcessedData;
    for j=1:length(data)
        sessType = [num2str(data{j}.Reward_L / data{j}.Reward_S) 'xFR' num2str(data{j}.NumPressRequired_S)];
        switch sessType
            case '2xFR6'
                sessInd = 1;
                FRind = 1;
            case '2xFR12'
                sessInd = 2;
                FRind = 2;
            case '5xFR6'
                sessInd = 3;
                FRind = 1;
            case '5xFR12'
                sessInd = 4;
                FRind = 2;
        end
        dayAbortedCount = 0;
        for k=1:data{j}.TotalTrialsCompleted
            count = count+1;
            mouseNames{count} = data{j}.AnimalName;
            relativeValue = (data{j}.Reward_L/data{j}.NumPressRequired_L(k)) - (data{j}.Reward_S/data{j}.NumPressRequired_S);
            if (isnan(data{j}.SideRewarded(k)))
                dayAbortedCount = dayAbortedCount + 1;
                if (data{j}.SideChosen(k) == 0)
                    abortPress = data{j}.LeverPressesByTrial(k); % - 1 to take out press that represents side choice
                    pressRequirement = data{j}.NumPressRequired_S; % - 1 to take out press that represents side choice
                    percentCompletedFR{sessInd} = [percentCompletedFR{sessInd} abortPress/pressRequirement];
                    
                    abort_countFR_all(FRind,abortPress) = abort_countFR_all(FRind,abortPress) + 1;
                    total_countFR_all(FRind,1:abortPress) = total_countFR_all(FRind,1:abortPress) + 1;
                    
                    abort_countFR{sessInd}(abortPress) = abort_countFR{sessInd}(abortPress) + 1;
                    total_countFR{sessInd}(1:abortPress) = total_countFR{sessInd}(1:abortPress) + 1;
                    
                    dataMatrix(:,count) = [data{j}.SideChosen(k); isnan(data{j}.SideRewarded(k)); data{j}.LeverPressesByTrial(k); abortPress/pressRequirement; relativeValue; sessInd; pressRequirement; k; j];
                elseif (data{j}.SideChosen(k) == 1)
                    abortPress = data{j}.LeverPressesByTrial(k);
                    pressRequirement = data{j}.NumPressRequired_L(k);
                    percentCompletedPR{sessInd} = [percentCompletedPR{sessInd} abortPress/pressRequirement];
                    
                    abort_countPR_all(pressRequirement,abortPress) = abort_countPR_all(pressRequirement,abortPress) + 1;
                    total_countPR_all(pressRequirement,1:abortPress) = total_countPR_all(pressRequirement,1:abortPress) + 1;
                    
                    abort_countPR{sessInd}(pressRequirement,abortPress) = abort_countPR{sessInd}(pressRequirement,abortPress) + 1;
                    total_countPR{sessInd}(pressRequirement,1:abortPress) = total_countPR{sessInd}(pressRequirement,1:abortPress) + 1;
                    
                    dataMatrix(:,count) = [data{j}.SideChosen(k); isnan(data{j}.SideRewarded(k)); data{j}.LeverPressesByTrial(k); abortPress/pressRequirement; relativeValue; sessInd; pressRequirement; k; j];
                elseif (isnan(data{j}.SideChosen(k)))
                    continue;
                else
                    error(['mouse ' mice{i} ' data{' num2str(j) '}.SideChosen(' num2str(k) ') is not valid'])
                end
            else
                if (data{j}.SideChosen(k) == 0)
                    pressRequirement = data{j}.NumPressRequired_S;
                    
                    total_countFR_all(FRind,1:pressRequirement) = total_countFR_all(FRind,1:pressRequirement) + 1;
                    total_countFR{sessInd} = total_countFR{sessInd} + 1;
                    
                    dataMatrix(:,count) = [data{j}.SideChosen(k); isnan(data{j}.SideRewarded(k)); data{j}.LeverPressesByTrial(k); 1; relativeValue; sessInd; pressRequirement; k; j];
                elseif (data{j}.SideChosen(k) == 1)
                    pressRequirement = data{j}.NumPressRequired_L(k);
                    
                    total_countPR_all(pressRequirement,1:pressRequirement) = total_countPR_all(pressRequirement,1:pressRequirement) + 1;
                    total_countPR{sessInd}(pressRequirement,1:pressRequirement) = total_countPR{sessInd}(pressRequirement,1:pressRequirement) + 1;
                    
                    dataMatrix(:,count) = [data{j}.SideChosen(k); isnan(data{j}.SideRewarded(k)); data{j}.LeverPressesByTrial(k); 1; relativeValue; sessInd; pressRequirement; k; j];
                elseif (isnan(data{j}.SideChosen(k)))
                    continue;
                else
                    error(['mouse ' mice{i} ' data{' num2str(j) '}.SideChosen(' num2str(k) ') is not valid'])
                end
            end
        end
        nAborted{sessInd} = [nAborted{sessInd} dayAbortedCount];
    end
end

dataTable = array2table(dataMatrix','VariableNames',{'SideChosen','WasAborted','LeverPressesByTrial','FracCompleted','RelativeValue','SessionIndex','PressRequirment','TrialNum','DayIndex'});
dataTable.mouseNames = mouseNames';

Pabort_PR_all = abort_countPR_all./total_countPR_all;
Pabort_FR_all = abort_countFR_all./total_countFR_all;
for i=1:4
    Pabort_PR{i} = abort_countPR{i}./total_countPR{i};
    Pabort_FR{i} = abort_countFR{i}./total_countFR{i};
end

save([saveDir 'percentCompletedFR.mat'],'percentCompletedFR','-mat')
save([saveDir 'percentCompletedPR.mat'],'percentCompletedPR','-mat')
save([saveDir 'Pabort_FR.mat'],'Pabort_FR','-mat')
save([saveDir 'Pabort_PR.mat'],'Pabort_PR','-mat')
save([saveDir 'Pabort_PR_all.mat'],'Pabort_PR_all','-mat')
save([saveDir 'Pabort_FR_all.mat'],'Pabort_FR_all','-mat')
save([saveDir 'dataMatrix.mat'],'dataMatrix','-mat')

%% Sigmoidal fit to side choice - 2 parameters
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
f = @(x,A,B) 1./(1 + exp(-A*x)) + B;
ft = fittype('1/(1 + exp(-A*x)) + B');
opts = fitoptions('METHOD','NonlinearLeastSquares');
F = fit(dataTable.RelativeValue,dataTable.SideChosen,ft,opts);
x = -1:.001:5;
plot(x,f(x,F.A,F.B))
xlabel('Relative Value (PR - FR)')
ylabel('P_{choose}(PR)')
figName='SideChoice_vs_relativeValue_sigmoidalFit';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')


figure;
As = zeros(1,4);
Bs = zeros(1,4);
bins = [-.4:.1:4.9; -.3:.1:5];
meanbins = mean(bins,1);
for i=1:4
    sessInds = find(dataTable.SessionIndex == i);
    for j=1:size(bins,2)
        binInds = find(dataTable.RelativeValue > bins(1,j) & dataTable.RelativeValue <= bins(2,j));
        curInds = intersect(sessInds,binInds);
        fracPRchosen{i}(j) = mean(dataTable.SideChosen(curInds));
    end
    F = fit(dataTable.RelativeValue(sessInds),dataTable.SideChosen(sessInds),ft,opts);
    As(i) = F.A;
    Bs(i) = F.B;
    plot(x,f(x,F.A,F.B))
    hold on;
end
Fall = fit(dataTable.RelativeValue,dataTable.SideChosen,ft,opts);
plot(x,f(x,Fall.A,Fall.B))
legend({['2xFR6: \tau = ' num2str(As(1)) ' , b = ' num2str(Bs(1))],...
        ['2xFR12: \tau = ' num2str(As(2)) ' , b = ' num2str(Bs(2))],...
        ['5xFR6: \tau = ' num2str(As(3)) ' , b = ' num2str(Bs(3))],...
        ['5xFR12: \tau = ' num2str(As(4)) ' , b = ' num2str(Bs(4))],...
        ['All data: \tau = ' num2str(Fall.A) ' , b = ' num2str(F.B)]},'Location','southeast')
xlabel('Relative Value (PR - FR)')
ylabel('P_{choose}(PR)')
set(gcf,'Position',[10 10 1400 1200])
figName = 'SideChoice_vs_relativeValue_sigmoidalFit_by_SessionType';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

figure;
for i=1:4
    subplot(2,2,i)
    scatter(meanbins,fracPRchosen{i},20,'.')
    hold on;
    plot(x,f(x,As(i),Bs(i)))
    title(sessionTypes{i})
    ylabel('P_{choose}(PR)')
    xlabel('Relative Value (PR - FR)')
end
set(gcf,'Position',[10 10 1400 1200])
figName = 'SideChoice_vs_relativeValue_sigmoidalFit_by_SessionType_w_data';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

% zoomed in fits
clear fracPRchosen As Bs
figure;
As = zeros(1,4);
Bs = zeros(1,4);
bins = [-.4:.05:.45; -.35:.05:.5];
x = -.4:.001:.5;
meanbins = mean(bins,1);
zoomInds = find(dataTable.RelativeValue >= bins(1,1) & dataTable.RelativeValue <= bins(2,end));
for i=1:4
    sessInds = find(dataTable.SessionIndex == i);
    for j=1:size(bins,2)
        binInds = find(dataTable.RelativeValue > bins(1,j) & dataTable.RelativeValue <= bins(2,j));
        curInds = intersect(sessInds,binInds);
        fracPRchosen{i}(j) = mean(dataTable.SideChosen(curInds));
    end
    F = fit(dataTable.RelativeValue(intersect(sessInds,zoomInds)),dataTable.SideChosen(intersect(sessInds,zoomInds)),ft,opts);
    As(i) = F.A;
    Bs(i) = F.B;
    plot(x,f(x,F.A,F.B))
    hold on;
end
Fall = fit(dataTable.RelativeValue(zoomInds),dataTable.SideChosen(zoomInds),ft,opts);
plot(x,f(x,Fall.A,Fall.B))
legend({['2xFR6: \tau = ' num2str(As(1)) ' , b = ' num2str(Bs(1))],...
        ['2xFR12: \tau = ' num2str(As(2)) ' , b = ' num2str(Bs(2))],...
        ['5xFR6: \tau = ' num2str(As(3)) ' , b = ' num2str(Bs(3))],...
        ['5xFR12: \tau = ' num2str(As(4)) ' , b = ' num2str(Bs(4))],...
        ['All data: \tau = ' num2str(Fall.A) ' , b = ' num2str(F.B)]},'Location','southeast')
xlabel('Relative Value (PR - FR)')
ylabel('P_{choose}(PR)')
set(gcf,'Position',[10 10 1400 1200])
figName = 'SideChoice_vs_relativeValue_sigmoidalFit_by_SessionType_zoom';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

figure;
for i=1:4
    subplot(2,2,i)
    scatter(meanbins,fracPRchosen{i},20,'.')
    hold on;
    plot(x,f(x,As(i),Bs(i)))
    title(sessionTypes{i})
    xlabel('Relative Value (PR - FR)')
    ylabel('P_{choose}(PR)')
end
set(gcf,'Position',[10 10 1400 1200])
figName = 'SideChoice_vs_relativeValue_sigmoidalFit_by_SessionType_w_data_zoom';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%% Sigmoidal fit to side choice - 4 parameters
clear fracPRchosen As Bs
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
f = @(x,A,B,C,D) B + (C - B)./(1 + exp(-A*(x - D)));
ft = fittype('B + (C - B)./(1 + exp(-A*(x - D)))');
opts = fitoptions('METHOD','NonlinearLeastSquares');
F = fit(dataTable.RelativeValue,dataTable.SideChosen,ft,opts);
x = -1:.001:5;
plot(x,f(x,F.A,F.B,F.C,F.D))
xlabel('Relative Value (PR - FR)')
ylabel('P_{choose}(PR)')
figName='SideChoice_vs_relativeValue_sigmoidalFit_4params';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')


figure;
As = zeros(1,4);
Bs = zeros(1,4);
Cs = zeros(1,4);
Ds = zeros(1,4);
bins = [-.4:.1:4.9; -.3:.1:5];
meanbins = mean(bins,1);
for i=1:4
    sessInds = find(dataTable.SessionIndex == i);
    for j=1:size(bins,2)
        binInds = find(dataTable.RelativeValue > bins(1,j) & dataTable.RelativeValue <= bins(2,j));
        curInds = intersect(sessInds,binInds);
        fracPRchosen{i}(j) = mean(dataTable.SideChosen(curInds));
    end
    F = fit(dataTable.RelativeValue(sessInds),dataTable.SideChosen(sessInds),ft,opts);
    As(i) = F.A;
    Bs(i) = F.B;
    Cs(i) = F.C;
    Ds(i) = F.D;
    plot(x,f(x,F.A,F.B,F.C,F.D))
    hold on;
end
Fall = fit(dataTable.RelativeValue,dataTable.SideChosen,ft,opts);
plot(x,f(x,Fall.A,Fall.B,Fall.C,Fall.D))
legend({['2xFR6: \tau = ' num2str(As(1)) ' , y_{min} = ' num2str(Bs(1)) ' , y_{max} = ' num2str(Cs(1)) ' , x0 = ' num2str(Ds(1))],...
        ['2xFR12: \tau = ' num2str(As(2)) ' , y_{min} = ' num2str(Bs(2)) ' , y_{max} = ' num2str(Cs(2)) ' , x0 = ' num2str(Ds(2))],...
        ['5xFR6: \tau = ' num2str(As(3)) ' , y_{min} = ' num2str(Bs(3))  ' , y_{max} = ' num2str(Cs(3)) ' , x0 = ' num2str(Ds(3))],...
        ['5xFR12: \tau = ' num2str(As(4)) ' , y_{min} = ' num2str(Bs(4))  ' , y_{max} = ' num2str(Cs(4)) ' , x0 = ' num2str(Ds(4))],...
        ['All data: \tau = ' num2str(Fall.A) ' , y_{min} = ' num2str(Fall.B) ' , y_{max} = ' num2str(Fall.C) ' , x0 = ' num2str(Fall.D)]},...
        'Location','southeast')
xlabel('Relative Value (PR - FR)')
ylabel('P_{choose}(PR)')
set(gcf,'Position',[10 10 1400 1200])
figName = 'SideChoice_vs_relativeValue_sigmoidalFit_by_SessionType_4params';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

figure;
for i=1:4
    subplot(2,2,i)
    scatter(meanbins,fracPRchosen{i},20,'.')
    hold on;
    plot(x,f(x,As(i),Bs(i),Cs(i),Ds(i)))
    title(sessionTypes{i})
    ylabel('P_{choose}(PR)')
    xlabel('Relative Value (PR - FR)')
end
set(gcf,'Position',[10 10 1400 1200])
figName = 'SideChoice_vs_relativeValue_sigmoidalFit_by_SessionType_w_data_4params';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

% zoomed in fits
clear fracPRchosen As Bs
figure;
As = zeros(1,4);
Bs = zeros(1,4);
Cs = zeros(1,4);
Ds = zeros(1,4);
bins = [-.4:.05:.45; -.35:.05:.5];
x = -.4:.001:.5;
meanbins = mean(bins,1);
zoomInds = find(dataTable.RelativeValue >= bins(1,1) & dataTable.RelativeValue <= bins(2,end));
for i=1:4
    sessInds = find(dataTable.SessionIndex == i);
    for j=1:size(bins,2)
        binInds = find(dataTable.RelativeValue > bins(1,j) & dataTable.RelativeValue <= bins(2,j));
        curInds = intersect(sessInds,binInds);
        fracPRchosen{i}(j) = mean(dataTable.SideChosen(curInds));
    end
    F = fit(dataTable.RelativeValue(intersect(sessInds,zoomInds)),dataTable.SideChosen(intersect(sessInds,zoomInds)),ft,opts);
    As(i) = F.A;
    Bs(i) = F.B;
    Cs(i) = F.C;
    Ds(i) = F.D;
    plot(x,f(x,F.A,F.B,F.C,F.D))
    hold on;
end
Fall = fit(dataTable.RelativeValue(zoomInds),dataTable.SideChosen(zoomInds),ft,opts);
plot(x,f(x,Fall.A,Fall.B,Fall.C,Fall.D))
legend({['2xFR6: \tau = ' num2str(As(1)) ' , y_{min} = ' num2str(Bs(1)) ' , y_{max} = ' num2str(Cs(1)) ' , x0 = ' num2str(Ds(1))],...
        ['2xFR12: \tau = ' num2str(As(2)) ' , y_{min} = ' num2str(Bs(2)) ' , y_{max} = ' num2str(Cs(2)) ' , x0 = ' num2str(Ds(2))],...
        ['5xFR6: \tau = ' num2str(As(3)) ' , y_{min} = ' num2str(Bs(3))  ' , y_{max} = ' num2str(Cs(3)) ' , x0 = ' num2str(Ds(3))],...
        ['5xFR12: \tau = ' num2str(As(4)) ' , y_{min} = ' num2str(Bs(4))  ' , y_{max} = ' num2str(Cs(4)) ' , x0 = ' num2str(Ds(4))],...
        ['All data: \tau = ' num2str(Fall.A) ' , y_{min} = ' num2str(Fall.B) ' , y_{max} = ' num2str(Fall.C) ' , x0 = ' num2str(Fall.D)]},...
        'Location','southeast')
xlabel('Relative Value (PR - FR)')
ylabel('P_{choose}(PR)')
set(gcf,'Position',[10 10 1400 1200])
figName = 'SideChoice_vs_relativeValue_sigmoidalFit_by_SessionType_zoom_4params';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

figure;
for i=1:4
    subplot(2,2,i)
    scatter(meanbins,fracPRchosen{i},20,'.')
    hold on;
    plot(x,f(x,As(i),Bs(i),Cs(i),Ds(i)))
    title(sessionTypes{i})
    ylabel('P_{choose}(PR)')
    xlabel('Relative Value (PR - FR)')
end
set(gcf,'Position',[10 10 1400 1200])
figName = 'SideChoice_vs_relativeValue_sigmoidalFit_by_SessionType_w_data_zoom_4params';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
figure;
for i=1:4
    subplot(2,2,i)
    plot(Pabort_FR{i})
    xlabel('Lever press #')
    ylabel('P(abort)')
    title(sessionTypes{i})
end
figName = 'Pabort_FR_by_leverPress';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
plot(Pabort_PR_all(10:20,:)')
title('Requirement: 10-20')
xlabel('Lever Press #')
ylabel('P(abort)')
figName = 'Pabort_PR_all_10-20';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
plot(Pabort_PR_all(20:30,:)')
title('Requirement: 20-30')
xlabel('Lever Press #')
ylabel('P(abort)')
figName = 'Pabort_PR_all_20-30';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
plot(Pabort_PR_all(30:40,:)')
title('Requirement: 30-40')
xlabel('Lever Press #')
ylabel('P(abort)')
figName = 'Pabort_PR_all_30-40';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
for i=1:4
    subplot(2,2,i)
    histogram(percentCompletedPR{i})
    xlabel('% of requirement')
    ylabel('Probability density')
    title(sessionTypes{i})
end
figName = 'percentCompletedPR_by_sessionType';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
bins = [-.4:.1:4.9; -.3:.1:5];
nbins = length(bins);
meanbins = mean(bins,1);
PRinds = find(dataMatrix(1,:) == 1);
fracAborted = zeros(1,nbins);
ninds = zeros(1,nbins);
meanLeverPresses = zeros(1,nbins);
stdLeverPresses = zeros(1,nbins);
meanFracRequirement = zeros(1,nbins);
stdFracRequirement = zeros(1,nbins);
boxPlotYs = [];
boxPlotXs = [];
for i=1:size(bins,2)
    inds = find(dataMatrix(5,:) > bins(1,i) & dataMatrix(5,:) <= bins(2,i));
    prinds = intersect(PRinds,inds);
    fracAborted(i) = sum(dataMatrix(2,prinds))/length(prinds);
    ninds(i) = length(prinds);
    meanLeverPresses(i) = mean(dataMatrix(3,prinds)); 
    stdLeverPresses(i) = std(dataMatrix(3,prinds));
    meanFracRequirement(i) = mean(dataMatrix(4,prinds));
    stdFracRequirement(i) = std(dataMatrix(4,prinds));
    boxPlotYs = [boxPlotYs dataMatrix(4,prinds)];
    boxPlotXs = [boxPlotXs ones(1,length(prinds))*meanbins(i)];
end

%%
figure;
shadedErrorBar(meanbins,meanFracRequirement,stdFracRequirement./sqrt(ninds));
xlabel('Relative value (PR - FR)')
ylabel('Fraction of requirement before abort')
figName = 'fracRequirementCompleted_vs_relativeValue';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
notBoxPlot(boxPlotYs,boxPlotXs,'jitter',.05)
xlabel('Relative value (PR - FR)')
ylabel('Fraction of requirement before abort')
figName = 'fracRequirementCompleted_vs_relativeValue_boxplot';
xlim([min(bins(:)) max(bins(:))])
set(gcf,'Position',[10 10 1600 1000])
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'epfs')
saveas(gcf,[figDir figName],'png')

%% Fraction of trials aborted vs. relative value scatter + sigmoidal fit
PRtrials = dataTable.SideChosen == 1;
NaNinds = isnan(fracAborted);
ft = fittype('1 - 1/(1 + exp(-A*x)) + B');
F = fit(meanbins(~NaNinds)',fracAborted(~NaNinds)',ft);
Fall = fit(dataTable.RelativeValue(PRtrials),dataTable.WasAborted(PRtrials),ft);
f = @(x,A,B) 1 - 1./(1 + exp(-A*x)) + B;
x = -1:.001:5;
figure;
scatter(meanbins,fracAborted,'.'); 
hold on;
%plot(x,f(x,F.A,F.B))
plot(x,f(x,Fall.A,Fall.B))
legend({'Binned data',...
       ['logistic: \beta = ' num2str(Fall.A) ' , offset = ' num2str(Fall.B)]})
xlabel('Relative value (PR - FR)')
ylabel('Fraction of trials aborted')
figName = 'fracAborted_vs_relativeValue';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
print([figDir figName '.png'],'-dpng','-r600')

%%
figure;
shadedErrorBar(meanbins,meanLeverPresses,stdLeverPresses./sqrt(ninds))
xlabel('Relative value (PR - FR)')
ylabel('# of lever presses before abort')
figName = 'leverPresses_vs_relativeValue';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
bins = [-.4:.05:.45; -.35:.05:.5];
nbins = length(bins);
meanbins = mean(bins,1);
PRinds = find(dataMatrix(1,:) == 1);
FRinds = find(dataMatrix(1,:) == 0);
fracAborted = zeros(1,nbins);
fracAbortedFR = zeros(1,nbins);
ninds = zeros(1,nbins);
nindsFR = zeros(1,nbins);
meanLeverPresses = zeros(1,nbins);
stdLeverPresses = zeros(1,nbins);
meanLeverPressesFR = zeros(1,nbins);
stdLeverPressesFR = zeros(1,nbins);
meanFracRequirement = zeros(1,nbins);
stdFracRequirement = zeros(1,nbins);
meanFracRequirementFR = zeros(1,nbins);
stdFracRequirementFR = zeros(1,nbins);
boxPlotYs = [];
boxPlotXs = [];
wasAborted = [];
boxPlotYsFR = [];
boxPlotXsFR = [];
wasAbortedFR = [];
for i=1:size(bins,2)
    inds = find(dataMatrix(5,:) > bins(1,i) & dataMatrix(5,:) <= bins(2,i));
    prinds = intersect(PRinds,inds);
    frinds = intersect(FRinds,inds);
    fracAborted(i) = sum(dataMatrix(2,prinds))/length(prinds);
    fracAbortedFR(i) = sum(dataMatrix(2,frinds))/length(frinds);
    ninds(i) = length(prinds);
    nindsFR(i) = length(frinds);
    meanLeverPresses(i) = mean(dataMatrix(3,prinds));
    stdLeverPresses(i) = std(dataMatrix(3,prinds));
    meanLeverPressesFR(i) = mean(dataMatrix(3,frinds));
    stdLeverPressesFR(i) = std(dataMatrix(3,frinds));
    meanFracRequirement(i) = mean(dataMatrix(4,prinds));
    stdFracRequirement(i) = std(dataMatrix(4,prinds));
    meanFracRequirementFR(i) = mean(dataMatrix(4,frinds));
    stdFracRequirementFR(i) = std(dataMatrix(4,frinds));
    boxPlotYs = [boxPlotYs dataMatrix(4,prinds)];
    boxPlotXs = [boxPlotXs ones(1,length(prinds))*meanbins(i)];
    wasAborted = [wasAborted dataMatrix(2,prinds)];
    boxPlotYsFR = [boxPlotYsFR dataMatrix(4,frinds)];
    boxPlotXsFR = [boxPlotXsFR ones(1,length(frinds))*meanbins(i)];
    wasAbortedFR = [wasAbortedFR dataMatrix(2,frinds)];
end
P = polyfit(meanbins,fracAborted,1);
[rho,pval] = corr([meanbins' fracAborted']);

%%
figure;
shadedErrorBar(meanbins,meanFracRequirement,stdFracRequirement./sqrt(ninds));
xlabel('Relative value (PR - FR)')
ylabel('Fraction of requirement before abort')
figName = 'fracRequirementCompleted_vs_relativeValue_zoom';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
notBoxPlot(boxPlotYs,boxPlotXs,'jitter',.05)
xlabel('Relative value (PR - FR)')
ylabel('Fraction of requirement before abort')
xlim([min(bins(:)) max(bins(:))])
set(gcf,'Position',[10 10 1600 1000])
figName = 'fracRequirementCompleted_vs_relativeValue_zoom_boxplot';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
notBoxPlot(boxPlotYs(logical(wasAborted)),boxPlotXs(logical(wasAborted)),'jitter',.03)
xlabel('Relative value (PR - FR)')
ylabel('Fraction of requirement before abort')
title('Aborted trials only')
xlim([min(bins(:)) max(bins(:))])
set(gcf,'Position',[10 10 1600 1000])
figName = 'fracRequirementCompleted_vs_relativeValue_zoom_abortOnly_boxplot';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
notBoxPlot(boxPlotYsFR(logical(wasAbortedFR)),boxPlotXsFR(logical(wasAbortedFR)),'jitter',.03)
xlabel('Relative value (PR - FR)')
ylabel('Fraction of requirement before abort')
title('FR: Aborted trials only')
xlim([min(bins(:)) max(bins(:))])
set(gcf,'Position',[10 10 1600 1000])
figName = 'fracRequirementCompletedFR_vs_relativeValue_zoom_abortOnly_boxplot';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
scatter(meanbins,fracAborted,'.');
hold on;
plot(meanbins,meanbins*P(1) + P(2))
text(.2,.7,{['y = ' num2str(P(1)) 'x + ' num2str(P(2))], ...
    '\rho = ' num2str(rho(1,2))})
xlabel('Relative value (PR - FR)')
ylabel('Fraction of trials aborted')
figName = 'fracAborted_vs_relativeValue_zoom';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
print([figDir figName '.png'],'-dpng','-r600')

%%
figure;
shadedErrorBar(meanbins,meanLeverPresses,stdLeverPresses./sqrt(ninds))
xlabel('Relative value (PR - FR)')
ylabel('# of lever presses before abort')
figName = 'leverPresses_vs_relativeValue_zoom';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
bins = [-.4:.05:.45; -.35:.05:.5];
nbins = length(bins);
meanbins = mean(bins,1);
PRinds = find(dataMatrix(1,:) == 1);
FRinds = find(dataMatrix(1,:) == 0);
boxPlotYsPR = cell(1,4);
boxPlotXsPR = cell(1,4);
wasAbortedPR = cell(1,4);
boxPlotYsFR = cell(1,4);
boxPlotXsFR = cell(1,4);
wasAbortedFR = cell(1,4);
sessTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
for i=1:4
    for j=1:size(bins,2)
        sessInds = find(dataMatrix(6,:) == i);
        inds = find(dataMatrix(5,:) > bins(1,j) & dataMatrix(5,:) <= bins(2,j));
        prinds = intersect(sessInds,intersect(PRinds,inds));
        frinds = intersect(sessInds,intersect(FRinds,inds));
        boxPlotYsPR{i} = [boxPlotYsPR{i} dataMatrix(4,prinds)];
        boxPlotXsPR{i} = [boxPlotXsPR{i} ones(1,length(prinds))*meanbins(j)];
        wasAbortedPR{i} = [wasAbortedPR{i} dataMatrix(2,prinds)];
        boxPlotYsFR{i} = [boxPlotYsFR{i} dataMatrix(4,frinds)];
        boxPlotXsFR{i} = [boxPlotXsFR{i} ones(1,length(frinds))*meanbins(j)];
        wasAbortedFR{i} = [wasAbortedFR{i} dataMatrix(2,frinds)];
    end
end

%%
figure;
for i=1:4
    subplot(2,2,i)
    notBoxPlot(boxPlotYsPR{i},boxPlotXsPR{i},'jitter',.04)
    xlabel('Relative value (PR - FR)')
    ylabel('Fraction of requirement before abort')
    xlim([min(bins(:)) max(bins(:))])
end
suptitle('PR trials: all trials')
set(gcf,'Position',[10 10 1600 1400])
figName = 'fracRequirementCompleted_vs_relativeValue_bySessionType_zoom';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
for i=1:4
    subplot(2,2,i)
    notBoxPlot(boxPlotYsPR{i}(logical(wasAbortedPR{i})),boxPlotXsPR{i}(logical(wasAbortedPR{i})),'jitter',.04)
    xlabel('Relative value (PR - FR)')
    ylabel('Fraction of requirement before abort')
    xlim([min(bins(:)) max(bins(:))])
end
suptitle('PR trials: aborted trials only')
set(gcf,'Position',[10 10 1600 1400])
figName = 'fracRequirementCompleted_vs_relativeValue_bySessionType_zoom_abortOnly';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
for i=1:4
    subplot(2,2,i)
    notBoxPlot(boxPlotYsFR{i},boxPlotXsFR{i},'jitter',.04)
    xlabel('Relative value (PR - FR)')
    ylabel('Fraction of requirement before abort')
    xlim([min(bins(:)) max(bins(:))])
end
suptitle('FR trials: all trials')
set(gcf,'Position',[10 10 1600 1400])
figName = 'fracRequirmentCompletedFR_vs_relativeValue_bySessionType_zoom';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%
figure;
for i=1:4
    subplot(2,2,i)
    notBoxPlot(boxPlotYsFR{i}(logical(wasAbortedFR{i})),boxPlotXsFR{i}(logical(wasAbortedFR{i})),'jitter',.04)
    xlabel('Relative value (PR - FR)')
    ylabel('Fraction of requirement before abort')
    xlim([min(bins(:)) max(bins(:))])
end
suptitle('FR trials: aborted trials only')
set(gcf,'Position',[10 10 1600 1400])
figName = 'fracRequirementCompletedFR_vs_relativeValue_bySessionType_zoom_abortOnly';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%% Hazard functions
% FR first
clear percentCompletedCDFs_FR percentCompletedCDFs_FR_allTrials pdfs_FR pdfs_FR_allTrials survivor_FR survivor_FR_allTrials abortConditionalProb_FR abortConditionalProb_FR_allTrials
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
for i=1:4
    if (i == 1 || i == 3)
        req = 6;
    else
        req = 12;
    end
    sessInds = find(dataTable.SessionIndex == i);
    frInds = find(dataTable.SideChosen == 0);
    curInds = intersect(sessInds,frInds);
    [n,edges] = histcounts(percentCompletedFR{i},(1:(req+1))/req);
    [n_allTrials,~] = histcounts(dataTable.FracCompleted(curInds),(1:(req+1))/req);
    n_allTrials(end) = 0;
    %n_allTrials = n_allTrials./length(curInds);
    percentCompletedCDFs_FR{i} = cumsum(n)/sum(n);
    percentCompletedCDFs_FR_allTrials{i} = cumsum(n_allTrials)/length(curInds); %sum(n_allTrials);
    pdfs_FR{i} = n./sum(n);
    pdfs_FR_allTrials{i} = n_allTrials./length(curInds);
    survivor_FR{i} = 1 - percentCompletedCDFs_FR{i};
    survivor_FR_allTrials{i} = 1 - percentCompletedCDFs_FR_allTrials{i};
end

for i=1:4
    for j=1:(length(survivor_FR{i}) - 1)
        abortConditionalProb_FR{i}(j) = (survivor_FR{i}(j) - survivor_FR{i}(j+1))/survivor_FR{i}(j);
    end
    for j=1:(length(survivor_FR_allTrials{i}) - 1)
        abortConditionalProb_FR_allTrials{i}(j) = (survivor_FR_allTrials{i}(j) - survivor_FR_allTrials{i}(j+1))/survivor_FR_allTrials{i}(j);
    end
end

figure;
for i=1:4
    subplot(2,2,i)
    plot(survivor_FR{i})
    ylabel('S(x)')
    xlabel('Lever press')
    title(sessionTypes{i})
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('FR trials')
figName = 'FR_trial_survival_functions';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

figure;
for i=1:4
    subplot(2,2,i)
    plot(pdfs_FR{i}./survivor_FR{i})
    xlabel('Lever press')
    ylabel('Hazard rate: p(x) / S(x)')
    if (i == 1 || i == 3)
        set(gca,'xtick',1:6)
    else
        set(gca,'xtick',1:12)
    end
    title(sessionTypes{i})
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('FR trials')
figName = 'FR_trial_hazard_rates';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

% Aborted FR trials only
figure;
for i=1:4
    subplot(2,2,i)
    plot(abortConditionalProb_FR{i})
    xlabel('Lever press')
    ylabel({'Conditional probability', 'of abort after press x'})
    ylim([0 1])
    if (i == 1 || i == 3)
        set(gca,'xtick',1:5)
    else
        set(gca,'xtick',1:11)
    end
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('Aborted FR trials')
figName = 'FR_trial_conditional_probability_of_abort_after_press_x';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

% All FR trials
figure;
for i=1:4
    subplot(2,2,i)
    plot(abortConditionalProb_FR_allTrials{i})
    xlabel('Lever press')
    ylabel({'Conditional probability', 'of abort after press x'})
    %ylim([0 1])
    if (i == 1 || i == 3)
        set(gca,'xtick',1:5)
    else
        set(gca,'xtick',1:11)
    end
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('All FR trials')
figName = 'FR_trial_conditional_probability_of_abort_after_press_x_allTrials';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PR trials
clear percentCompletedCDFs_PR percentCompletedCDFs_PR_allTrials pdfs_PR pdfs_PR_allTrials survivor_PR survivor_PR_allTrials abortConditionalProb_PR abortConditionalProb_PR_allTrials
clear percentCompletedCDFs_PR_allLowTrials percentCompletedCDFs_PR_allHighTrials
clear pdfs_PR_allLowTrials pdfs_PR_allHighTrials
clear survivor_PR_allLowTrials survivor_PR_allHighTrials
for i=1:4
    sessInds = find(dataTable.SessionIndex == i);
    prInds = find(dataTable.SideChosen == 1);
    lowValInds = find(dataTable.RelativeValue < -.05);
    highValInds = find(dataTable.RelativeValue > .2);
    curInds = intersect(sessInds,prInds);
    curLowInds = intersect(curInds,lowValInds);
    curHighInds = intersect(curInds,highValInds);
    [n,edges]=histcounts(percentCompletedPR{i},100);
    [n_allTrials,~] = histcounts(dataTable.FracCompleted(curInds),100);
    [n_allLowTrials,~] = histcounts(dataTable.FracCompleted(curLowInds),100);
    [n_allHighTrials,~] = histcounts(dataTable.FracCompleted(curHighInds),100);
    n_allTrials(end) = 0;
    n_allLowTrials(end) = 0;
    n_allHighTrials(end) = 0;
    percentCompletedCDFs_PR{i} = cumsum(n)/sum(n);
    percentCompletedCDFs_PR_allTrials{i} = cumsum(n_allTrials)/length(curInds);
    percentCompletedCDFs_PR_allLowTrials{i} = cumsum(n_allLowTrials)/length(curLowInds);
    percentCompletedCDFs_PR_allHighTrials{i} = cumsum(n_allHighTrials)/length(curHighInds);
    pdfs_PR{i} = n./sum(n);
    pdfs_PR_allTrials{i} = n_allTrials./length(curInds);
    pdfs_PR_allLowTrials{i} = n_allLowTrials./length(curLowInds);
    pdfs_PR_allHighTrials{i} = n_allHighTrials./length(curHighInds);
    survivor_PR{i} = 1 - percentCompletedCDFs_PR{i};
    survivor_PR_allTrials{i} = 1 - percentCompletedCDFs_PR_allTrials{i};
    survivor_PR_allLowTrials{i} = 1 - percentCompletedCDFs_PR_allLowTrials{i};
    survivor_PR_allHighTrials{i} = 1 - percentCompletedCDFs_PR_allHighTrials{i};
end

for i=1:4
    for j=1:(length(survivor_PR{i}) - 1)
        abortConditionalProb_PR{i}(j) = (survivor_PR{i}(j) - survivor_PR{i}(j+1))/survivor_PR{i}(j);
    end
    for j=1:(length(survivor_PR_allTrials{i}) - 1)
        abortConditionalProb_PR_allTrials{i}(j) = (survivor_PR_allTrials{i}(j) - survivor_PR_allTrials{i}(j+1))/survivor_PR_allTrials{i}(j);
    end
    for j=1:(length(survivor_PR_allLowTrials{i}) - 1)
        abortConditionalProb_PR_allLowTrials{i}(j) = (survivor_PR_allLowTrials{i}(j) - survivor_PR_allLowTrials{i}(j+1))/survivor_PR_allLowTrials{i}(j);
    end
    for j=1:(length(survivor_PR_allHighTrials{i}) - 1)
        abortConditionalProb_PR_allHighTrials{i}(j) = (survivor_PR_allHighTrials{i}(j) - survivor_PR_allHighTrials{i}(j+1))/survivor_PR_allHighTrials{i}(j);
    end
end

figure;
for i=1:4
    subplot(2,2,i)
    plot(survivor_PR{i})
    ylabel('S(x)')
    xlabel('% of trial completed')
    title(sessionTypes{i})
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('PR trials')
figName = 'PR_trial_survival_functions';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

figure;
for i=1:4
    subplot(2,2,i)
    plot(pdfs_PR{i}./survivor_PR{i})
    xlabel('% of trial completed')
    ylabel('Hazard rate: p(x) / S(x)')
    title(sessionTypes{i})
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('PR trials')
figName = 'PR_trial_hazard_rates';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

% Aborted Trials only
figure;
for i=1:4
    subplot(2,2,i)
    plot(abortConditionalProb_PR{i})
    xlabel('% of trial completed')
    ylabel({'Conditional probability of','abort after x% of trial'})
    ylim([0 1])
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('Aborted PR trials only')
set(gcf,'Position',[10 10 1400 1200])
figName = 'PR_trial_conditional_probability_of_abort_after_press_x';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

% All trials
figure;
for i=1:4
    subplot(2,2,i)
    plot(abortConditionalProb_PR_allTrials{i})
    xlabel('% of trial completed')
    ylabel({'Conditional probability of','abort after x% of trial'})
    %ylim([0 1])
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('All PR trials')
set(gcf,'Position',[10 10 1400 1200])
figName = 'PR_trial_conditional_probability_of_abort_after_press_x_allTrials';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

% All trials High vs. Low value difference
figure;
for i=1:4
    subplot(2,2,i)
    plot(abortConditionalProb_PR_allLowTrials{i})
    hold on;
    plot(abortConditionalProb_PR_allHighTrials{i})
    if (i == 1)
        legend({'Q_{PR} < Q_{FR}','Q_{PR} > Q_{FR}'})
    end
    xlabel('% of trial completed')
    ylabel({'Conditional probability of','abort after x% of trial'})
    %ylim([0 1])
    title(sessionTypes{i})
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('All PR trials. High vs. low value')
set(gcf,'Position',[10 10 1400 1200])
figName = 'PR_trial_conditional_probability_of_abort_after_press_x_allTrials_High_vs_Low_value';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

% PR trials less bins
clear abortConditionalProb_PR
for i=1:4
    [n,edges]=histcounts(percentCompletedPR{i},20);
    percentCompletedCDFs_PR{i} = cumsum(n)/sum(n);
    pdfs_PR{i} = n./sum(n);
    survivor_PR{i} = 1 - percentCompletedCDFs_PR{i};
end

for i=1:4
    for j=1:(length(survivor_PR{i}) - 1)
        abortConditionalProb_PR{i}(j) = (survivor_PR{i}(j) - survivor_PR{i}(j+1))/survivor_PR{i}(j);
    end
end

figure;
for i=1:4
    subplot(2,2,i)
    plot(1 - percentCompletedCDFs_PR{i})
    ylabel('S(x)')
    xlabel('% of trial completed')
    title(sessionTypes{i})
    set(gca,'xtick',0:20,'xticklabels',0:5:100)
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('PR trials')
figName = 'PR_trial_survival_functions_20bins';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

maxHazard = 0;
for i=1:4
    hazardRates{i} = pdfs_PR{i}./(1 - percentCompletedCDFs_PR{i});
    hazardRates{i}(isinf(hazardRates{i})) = NaN;
    if (max(hazardRates{i}) > maxHazard)
        maxHazard = max(hazardRates{i});
    end
end
figure;
for i=1:4
    subplot(2,2,i)
    plot(hazardRates{i})
    xlabel('% of trial completed')
    ylabel('Hazard rate: p(x) / S(x)')
    title(sessionTypes{i})
    set(gca,'xtick',0:20,'xticklabels',0:5:100)
    ylim([0 maxHazard])
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('PR trials')
figName = 'PR_trial_hazard_rates_20bins';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')

figure;
for i=1:4
    subplot(2,2,i)
    plot(abortConditionalProb_PR{i})
    xlabel('% of trial completed')
    ylabel({'Conditional probability of','abort after x% of trial'})
    ylim([0 1])
    set(gca,'xtick',0:20,'xticklabels',0:5:100)
end
set(gcf,'Position',[10 10 1400 1200])
suptitle('PR trials')
set(gcf,'Position',[10 10 1400 1200])
figName = 'PR_trial_conditional_probability_of_abort_after_press_x_20bins';
saveas(gcf,[figDir figName],'fig')
saveas(gcf,[figDir figName],'eps')
saveas(gcf,[figDir figName],'png')