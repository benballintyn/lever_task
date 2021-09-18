function [scores,allParams] = plotResults(datadir,top_k,varargin)
p=inputParser;
addRequired(p,'datadir',@ischar)
addRequired(p,'top_k',@isnumeric)
addParameter(p,'single_dir',false,@islogical)
addParameter(p,'fitStatsOnly',false,@islogical)
addParameter(p,'fitMethod','logprob_joint',@ischar)
parse(p,datadir,top_k,varargin{:})

NL_observed = load('optimality/WT/NL_observed.mat'); NL_observed=NL_observed.NL_observed;
nTrialsAborted = load('optimality/WT/nTrialsAborted.mat'); nTrialsAborted=nTrialsAborted.nTrialsAborted;
EoR_optimalities = load('optimality/WT/EoR_optimalities.mat'); EoR_optimalities=EoR_optimalities.EoR_optimalities;
percentCompletedPR_allTrials = load('driftRL/abort_analysis/percentCompletedPR_allTrials.mat'); percentCompletedPR_allTrials = percentCompletedPR_allTrials.percentCompletedPR_allTrials;

sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
if (p.Results.single_dir)
    score = load([datadir '/objective_score.mat']); score=score.score;
    numLR = load([datadir '/numLR.mat']); numLR=numLR.numLR;
    numAborted = load([datadir '/numAborted.mat']); numAborted=numAborted.numAborted;
    performanceEoR = load([datadir '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;
    percentCompletedPR = load([datadir '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
    params = load([datadir '/params.mat']); params=params.params;
    figure;
    for i=1:16
        subplot(4,4,i)
        if (i <= 4)
            histogram(NL_observed{i},'normalization','pdf','binwidth',2);
            hold on;
            histogram(numLR{i},'normalization','pdf','binwidth',2);
            xlabel('# of completed PR trials','FontSize',15,'FontWeight','bold')
            title([sessionTypes{i}])
        elseif (i > 4 & i <= 8)
            curInd = i - 4;
            histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2);
            hold on;
            histogram(numAborted{curInd},'normalization','pdf','binwidth',2);
            xlabel('# of trials aborted','FontSize',15,'FontWeight','bold')
        elseif (i > 8 & i <= 12)
            curInd = i - 8;
            histogram(EoR_optimalities{curInd},'normalization','pdf','binwidth',.02);
            hold on;
            histogram(performanceEoR{curInd},'normalization','pdf','binwidth',.02);
            xlabel('EoR optimality','FontSize',15,'FontWeight','bold')
        else
            curInd = i - 12;
            histogram(percentCompletedPR_allTrials{curInd},'normalization','pdf','binwidth',.02)
            hold on;
            histogram(percentCompletedPR{curInd},'normalization','pdf','binwidth',.02)
            xlabel('% of PR trial completed','FontSize',15,'FontWeight','bold')
        end
    end
    suptitle(num2str(params))
    set(gcf,'Position',[10 10 1400 1200])
    scores = score;
    allParams = params;
else
    dirContents = dir(datadir);
    dirContents = dirContents(~ismember({dirContents.name},{'..','.'}));
    dirContents = dirContents([dirContents.isdir]);
    dirContents
    for i=1:length(dirContents)
        folderNums(i) = str2num(dirContents(i).name);
    end
    [~,folderOrder] = sort(folderNums);
    dirContents = dirContents(folderOrder);

    for i=1:length(dirContents)
        s = load([datadir '/' dirContents(i).name '/objective_score.mat']); s=s.score;
        params = load([datadir '/' dirContents(i).name '/params.mat']); params = params.params;
        scores(i) = s;
        allParams(i,:) = params;
    end
    
    if (p.Results.fitStatsOnly)
        [~,inds] = sort(scores);
        for i=1:top_k
            numLR = load([datadir '/' dirContents(inds(i)).name '/numLR.mat']); numLR=numLR.numLR;
            numAborted = load([datadir '/' dirContents(inds(i)).name '/numAborted.mat']); numAborted=numAborted.numAborted;
            performanceEoR = load([datadir '/' dirContents(inds(i)).name '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;
            percentCompletedPR = load([datadir '/' dirContents(inds(i)).name '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
            params = load([datadir '/' dirContents(inds(i)).name '/params.mat']); params=params.params;
            figure;
            if (strcmp(p.Results.fitMethod,'logprob_joint'))
                for j=1:12
                    subplot(3,4,j)
                    if (j <= 4)
                        histogram(NL_observed{j},'normalization','pdf','binwidth',2);
                        hold on;
                        histogram(numLR{j},'normalization','pdf','binwidth',2);
                        xlabel('# of completed PR trials','FontSize',15,'FontWeight','bold')
                        title([sessionTypes{j}])
                    elseif (j > 4 & j <= 8)
                        curInd = j - 4;
                        histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2);
                        hold on;
                        histogram(numAborted{curInd},'normalization','pdf','binwidth',2);
                        xlabel('# of trials aborted','FontSize',15,'FontWeight','bold')
                    else
                        curInd = j - 8;
                        histogram(EoR_optimalities{curInd},'normalization','pdf','binwidth',.02);
                        hold on;
                        histogram(performanceEoR{curInd},'normalization','pdf','binwidth',.02);
                        xlabel('EoR optimality','FontSize',15,'FontWeight','bold')
                    end
                end
                curtitle = {num2str(params),['index = ' num2str(inds(i)) ' score = ' num2str(scores(inds(i)))]};
                suptitle(curtitle)
                set(gcf,'Position',[10 10 1400 1200])
            elseif (strcmp(p.Results.fitMethod,'logprob_independent'))
                for j=1:16
                    subplot(4,4,j)
                    if (j <= 4)
                        histogram(NL_observed{j},'normalization','pdf','binwidth',2);
                        hold on;
                        histogram(numLR{j},'normalization','pdf','binwidth',2);
                        xlabel('# of completed PR trials','FontSize',15,'FontWeight','bold')
                        title([sessionTypes{j}])
                    elseif (j > 4 & j <= 8)
                        curInd = j - 4;
                        histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2);
                        hold on;
                        histogram(numAborted{curInd},'normalization','pdf','binwidth',2);
                        xlabel('# of trials aborted','FontSize',15,'FontWeight','bold')
                    elseif (j > 8 & j <= 12)
                        curInd = j - 8;
                        histogram(EoR_optimalities{curInd},'normalization','pdf','binwidth',.02);
                        hold on;
                        histogram(performanceEoR{curInd},'normalization','pdf','binwidth',.02);
                        xlabel('EoR optimality','FontSize',15,'FontWeight','bold')
                    else
                        curInd = j - 12;
                        histogram(percentCompletedPR_allTrials{curInd},'normalization','pdf','binwidth',.02)
                        hold on;
                        histogram(percentCompletedPR{curInd},'normalization','pdf','binwidth',.02)
                        xlabel('% of PR trial completed','FontSize',15,'FontWeight','bold')
                    end
                end
                curtitle = {num2str(params),['index = ' num2str(inds(i)) ' score = ' num2str(scores(inds(i)))]};
                suptitle(curtitle)
                set(gcf,'Position',[10 10 1400 1200])
            end
        end
    else
        [~,inds] = sort(scores);
        for i=1:top_k
            numLR = load([datadir '/' dirContents(inds(i)).name '/numLR.mat']); numLR=numLR.numLR;
            numAborted = load([datadir '/' dirContents(inds(i)).name '/numAborted.mat']); numAborted=numAborted.numAborted;
            performanceEoR = load([datadir '/' dirContents(inds(i)).name '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;
            percentCompletedPR = load([datadir '/' dirContents(inds(i)).name '/percentCompletedPR.mat']); percentCompletedPR=percentCompletedPR.percentCompletedPR;
            params = load([datadir '/' dirContents(inds(i)).name '/params.mat']); params=params.params;
            figure;
            for j=1:16
                subplot(4,4,j)
                if (j <= 4)
                    histogram(NL_observed{j},'normalization','pdf','binwidth',2);
                    hold on;
                    histogram(numLR{j},'normalization','pdf','binwidth',2);
                    xlabel('# of completed PR trials','FontSize',15,'FontWeight','bold')
                    title([sessionTypes{j}])
                elseif (j > 4 & j <= 8)
                    curInd = j - 4;
                    histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2);
                    hold on;
                    histogram(numAborted{curInd},'normalization','pdf','binwidth',2);
                    xlabel('# of trials aborted','FontSize',15,'FontWeight','bold')
                elseif (j > 8 & j <= 12)
                    curInd = j - 8;
                    histogram(EoR_optimalities{curInd},'normalization','pdf','binwidth',.02);
                    hold on;
                    histogram(performanceEoR{curInd},'normalization','pdf','binwidth',.02);
                    xlabel('EoR optimality','FontSize',15,'FontWeight','bold')
                else
                    curInd = j - 12;
                    histogram(percentCompletedPR_allTrials{curInd},'normalization','pdf','binwidth',.02)
                    hold on;
                    histogram(percentCompletedPR{curInd},'normalization','pdf','binwidth',.02)
                    xlabel('% of PR trial completed','FontSize',15,'FontWeight','bold')
                end
            end
            curtitle = {num2str(params),['index = ' num2str(inds(i)) ' score = ' num2str(scores(inds(i)))]};
            suptitle(curtitle)
            set(gcf,'Position',[10 10 1400 1200])
        end
    end
    
    if (~p.Results.fitStatsOnly)
        figure;
        plot(scores)

        figure;
        plot(allParams)
        
        [minScore,bestScoreInd] = min(scores);
        bestParams = allParams(bestScoreInd,:);
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(allParams);
        cmap = jet;
        scoreColorInds = round(((scores - minScore)/abs(minScore))*254 + 1);
        scoreColorInds(scoreColorInds > 256) = 256;
        figure;
        subplot(2,2,1)
        scatter3(allParams*COEFF(:,1),allParams*COEFF(:,2),allParams*COEFF(:,3),100,cmap(scoreColorInds,:),'.')
        hold on;
        plot3(bestParams*COEFF(:,1),bestParams*COEFF(:,2),bestParams*COEFF(:,3),'r*','MarkerSize',20)
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')

        subplot(2,2,2)
        plot(EXPLAINED)

        subplot(2,2,3)
        imagesc(COEFF)
    end
end
end

