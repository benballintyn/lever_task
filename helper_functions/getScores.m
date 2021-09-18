function [scores,allParams] = getScores(datadir,varargin)
p=inputParser;
addRequired(p,'datadir',@ischar)
addParameter(p,'useSave',false,@islogical)
parse(p,datadir,varargin{:})

dirContents = dir(datadir);
dirContents = dirContents(~ismember({dirContents.name},{'..','.'}));
dirContents = dirContents([dirContents.isdir]);
for i=1:length(dirContents)
    folderNums(i) = str2num(dirContents(i).name);
end
[~,folderOrder] = sort(folderNums);
dirContents = dirContents(folderOrder);

if (p.Results.useSave)
    if (exist([datadir '/allScores.mat'],'file'))
        scores = load([datadir '/allScores.mat']); scores=scores.scores;
        disp('Saved scores loaded')
    else
        disp('Saved scores not found. Loading...')
        for i=1:length(dirContents)
            s = load([datadir '/' dirContents(i).name '/objective_score.mat']); s=s.score;
            scores(i) = s;
        end
        save([datadir '/allScores.mat'],'scores','-mat')
        disp('Saved allScores')
    end
    if (exist([datadir '/allParams.mat']))
        allParams = load([datadir '/allParams.mat']); allParams=allParams.allParams;
        disp('Saved allParams loaded')
    else
        for i=1:length(dirContents)
            params = load([datadir '/' dirContents(i).name '/params.mat']); params = params.params;
            allParams(i,:) = params;
        end
        save([datadir '/allParams.mat'],'allParams','-mat')
        disp('Saved allParams')
    end
else
    for i=1:length(dirContents)
        s = load([datadir '/' dirContents(i).name '/objective_score.mat']); s=s.score;
        params = load([datadir '/' dirContents(i).name '/params.mat']); params = params.params;
        scores(i) = s;
        allParams(i,:) = params;
    end
    save([datadir '/allScores.mat'],'scores','-mat')
    save([datadir '/allParams.mat'],'allParams','-mat')
end
end

