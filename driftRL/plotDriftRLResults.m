function [] = plotDriftRLResults(datadir,top_k,varargin)
p=inputParser;
addRequired(p,'datadir',@ischar)
addRequired(p,'top_k',@isnumeric)
addParameter(p,'single_dir',false,@islogical)
parse(p,datadir,top_k,varargin{:})

NL_observed = load('optimality/WT/NL_observed.mat'); NL_observed=NL_observed.NL_observed;
nTrialsAborted = load('optimality/WT/nTrialsAborted.mat'); nTrialsAborted=nTrialsAborted.nTrialsAborted;


if (p.Results.single_dir)
    score = load([datadir '/objective_score.mat']); score=score.score;
    numLR = load([datadir '/numLR.mat']); numLR=numLR.numLR;
    numAborted = load([datadir '/numAborted.mat']); numAborted=numAborted.numAborted;
    params = load([datadir '/params.mat']); params=params.params;
    figure;
    for i=1:8
        subplot(2,4,i)
        if (i <= 4)
            histogram(NL_observed{i},'normalization','pdf','binwidth',2);
            hold on;
            histogram(numLR{i},'normalization','pdf','binwidth',2);
        else
            curInd = i - 4;
            histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2);
            hold on;
            histogram(numAborted{curInd},'normalization','pdf','binwidth',2);
        end
    end
    suptitle(num2str(params))
    set(gcf,'Position',[10 10 1400 1200])
    
else
    dirContents = dir(datadir);
    dirContents = dirContents(~ismember({dirContents.name},{'..','.'}));
    dirContents = dirContents([dirContents.isdir]);

    for i=1:length(dirContents)
        s = load([datadir '/' dirContents(i).name '/objective_score.mat']); s=s.score;
        scores(i) = s;
    end
    [~,inds] = sort(scores);

    for i=1:top_k
        numLR = load([datadir '/' dirContents(inds(i)).name '/numLR.mat']); numLR=numLR.numLR;
        numAborted = load([datadir '/' dirContents(inds(i)).name '/numAborted.mat']); numAborted=numAborted.numAborted;
        params = load([datadir '/' dirContents(inds(i)).name '/params.mat']); params=params.params;
        figure;
        for j=1:8
            subplot(2,4,j)
            if (j <= 4)
                histogram(NL_observed{j},'normalization','pdf','binwidth',2);
                hold on;
                histogram(numLR{j},'normalization','pdf','binwidth',2);
            else
                curInd = j - 4;
                histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2);
                hold on;
                histogram(numAborted{curInd},'normalization','pdf','binwidth',2);
            end
        end
        curtitle = {num2str(params),['index = ' num2str(inds(i))]};
        suptitle(curtitle)
        set(gcf,'Position',[10 10 1400 1200])
    end
end
end

