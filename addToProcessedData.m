% addToProcessedData
mice = {'AZ04','AZ08','HP01','HP02','HP03','HP04','MA01','NS07','NS09','NS10'};
for i=1:length(mice) % for each animal
    ProcessedData = load(['processed_data/' mice{i} '_ProcessedData.mat']); ProcessedData=ProcessedData.ProcessedData;
    filenames = dir(['raw_data/' mice{i} '*']);
    if (length(filenames) ~= length(ProcessedData))
        error('# of raw files does not equal length of ProcessedData')
    end
    for j=1:length(filenames) % for each day
        rawData = load(['raw_data/' filenames(j).name]); rawData=rawData.SessionData;
        for k=1:rawData.nTrials % for each trial
            curFieldNames = fieldnames(rawData.RawEvents.Trial{k}.States);
            if (sum(isnan(rawData.RawEvents.Trial{k}.States.(EndSession))) > 0)
                ProcessedData{j}.LeverPressesByTrial(k) = 0;
            else
                nPresses = 1;
                for l=1:length(curFieldNames) % for each event field
                    if (contains(curFieldNames{l},'Debounce'))
                        if (sum(isnan(rawData.RawEvents.Trial{k}.States.(curFieldNames{l}))) == 0)
                            nPresses = nPresses + 1;
                        end
                    end
                end
                ProcessedData{j}.LeverPressesByTrial(k) = nPresses;
            end
            disp(['Done with mouse ' mice{i} ' day ' num2str(j) ' trial ' num2str(k)])
        end
    end
    save(['processed_data/' mice{i} '_ReProcessedData.mat'],'ProcessedData','-mat')
end
    
