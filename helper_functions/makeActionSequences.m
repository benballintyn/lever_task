function [actionSeqs] = makeActionSequences(animal)
data = load(['processed_data/' animal '_ReProcessedData.mat']); data=data.ProcessedData;
nDays = length(data);
actionSeqs = cell(1,nDays);
for i=1:nDays
    actionSeqs{i} = [];
    for t=1:data{i}.TotalTrialsCompleted
        if (data{i}.SideChosen(t) == 0)
            actionSeqs{i} = [actionSeqs{i} 1];
        elseif (data{i}.SideChosen(t) == 1)
            actionSeqs{i} = [actionSeqs{i} 2];
        elseif (isnan(data{i}.SideChosen(t)))
            if (t == data{i}.TotalTrialsCompleted)
                continue;
            else
                error([animal ' day ' num2str(i) ' trial ' num2str(t) ' Side chosen is not valid'])
            end
        end
        actionSeqs{i} = [actionSeqs{i} ones(1,(data{i}.LeverPressesByTrial(t)-1))];
        if (isnan(data{i}.SideRewarded(t)))
            actionSeqs{i} = [actionSeqs{i} 2];
        end
    end
end

end

