clear all;
mice = {'AZ04','AZ08','HP01','HP02','HP03','HP04','MA01','NS07','NS09','NS10'};

count = 0;
for i=1:length(mice)
    d = load(['processed_data/' mice{i} '_ReProcessedData.mat']); d=d.ProcessedData;
    for j=1:length(d)
        for k=1:d{j}.TotalTrialsCompleted
            if (d{j}.SideChosen(k) == 1)
                for l=1:d{j}.LeverPressesByTrial(k)
                    count=count+1;
                    vals(count,1) = l;
                    vals(count,2) = d{j}.NumPressRequired_L(k) - l;
                    if (isnan(d{j}.SideRewarded(k)))
                        vals(count,3) = 0;
                    else
                        vals(count,3) = 1;
                    end
                end
            end
        end
    end
end
binMins = 0:5:50;
binMaxs = 5:5:55;
rem = 1:50;
cmap=jet;
cinds = floor(linspace(1,size(cmap,1),length(binMins)));
for i=1:length(binMins)
    for j=rem
        inds1 = find(vals(:,1) > binMins(i) & vals(:,1) < binMaxs(i));
        inds2 = find(vals(:,2)  == j);
        inds = intersect(inds1,inds2);
        if (~isempty(inds))
            rSum = sum(vals(inds,3))/length(inds);
            p(i,j) = rSum;
        end
    end
end
figure;
for i=1:size(p,1)
    scatter(rem,p(i,:),20,cmap(cinds(i),:),'filled')
    fits{i} = polyfit(rem,p(i,:),1);
    hold on;
end

figure;
for i=1:length(fits)
    plot(rem,rem*fits{i}(1) + fits{i}(2))
    hold on;
end