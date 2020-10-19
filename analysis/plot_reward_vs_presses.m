mice = {'AZ04','AZ08','HP01','HP02','HP03','HP04','MA01','NS07','NS09','NS10'};
count = 0;
cmap = jet;
figure;
for i=1:length(mice)
    data = load(['processed_data/WT/' mice{i} '_ReProcessedData.mat']);
    data=data.ProcessedData;
    
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
        plot(cumsum(data{j}.LeverPressesByTrial),'Color',cmap(ind*60,:)); hold on;
        count=count+1;
        totalR(count) = data{j}.TotalRewardCollected;
        totalPress(count) = sum(data{j}.LeverPressesByTrial);
        sesstype(count) = ind;
    end
end
figure;
scatter(totalPress,totalR,20,cmap(sesstype*60,:),'.')