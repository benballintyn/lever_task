clear all;
mice = {'AZ04','AZ08','HP01','HP02','HP03','HP04','MA01','NS07','NS09','NS10'};
count = 0;
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
        for k = 1:data{j}.TotalTrialsCompleted
            if (~isnan(data{j}.SideRewarded(k)))
                count = count + 1;
                %X(count,1) = SR;
                X(count,1) = LR;
                X(count,2) = Ps;
                X(count,3) = k;
                Y(count) = data{j}.SideChosen(k);
            end
        end
    end
end
[bLR,devLR,statsLR] = glmfit(X,Y','binomial');
[bSR,devSR,statsSR] = glmfit(X,(~Y)','binomial');
figure;
for i=1:4
    switch i
        case 1
            SR = 3; LR = 6; Ps = 6;
        case 2
            SR = 3; LR = 6; Ps = 12;
        case 3
            SR = 3; LR = 15; Ps = 6;
        case 4
            SR = 3; LR = 15; Ps = 12;
    end
    Xfit = [ones(800,1)*LR ones(800,1)*Ps (1:800)'];
    [ypredLR,yloLR,yhiLR] = glmval(bLR,Xfit,'logit',statsLR);
    [ypredSR,yloSR,yhiSR] = glmval(bSR,Xfit,'logit',statsSR);
    subplot(2,2,i)
    h1=shadedErrorBar(1:800,ypredLR',[yloLR'; yhiLR'],'lineprops','r');
    hold on;
    h2=shadedErrorBar(1:800,ypredSR',[yloSR'; yhiSR'],'lineprops','b');
    legend([h1.mainLine h2.mainLine],{'P(LR)','P(SR)'})
    xlabel('Trial #'); ylabel('P(side_t = LR/SR)')
    set(gca,'xtick',[1 200 400 600 800],'xticklabels',{'1','','','','800'},'ytick',[0 .2 .4 .6 .8 1],'yticklabels',{'0','','','','','1'})
end
