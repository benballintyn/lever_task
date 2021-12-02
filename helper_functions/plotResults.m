function [scores,allParams] = plotResults(datadir,top_k,varargin)
p=inputParser;
addRequired(p,'datadir',@ischar)
addRequired(p,'top_k',@isnumeric)
addParameter(p,'single_dir',false,@islogical)
addParameter(p,'fitStatsOnly',false,@islogical)
addParameter(p,'fitMethod','logprob_joint',@ischar)
addParameter(p,'Only120Trials',false,@islogical)
addParameter(p,'histograms_only',false,@islogical)
parse(p,datadir,top_k,varargin{:})

NUM_LR_VALS = 0:.01:100;
NUM_ABORTED_VALS = 0:1:60;
PERFORMANCE_EOR_VALS = 0:.01:1;
PERCENT_COMPLETED_PR_VALS = 0:.01:1;

NUM_LR_BANDWIDTH = 5;
NUM_ABORTED_BANDWIDTH = 5;
PERFORMANCE_EOR_BANDWIDTH = .05;
PERCENT_COMPLETED_PR_BANDWIDTH = .05;

if (p.Results.Only120Trials)
    mouseDataDir = '~/phd/lever_task/analysis_120_trials';
    NL_observed = load([mouseDataDir '/NL_observed.mat']); NL_observed=NL_observed.NL_observed;
    nTrialsAborted = load([mouseDataDir '/nTrialsAborted.mat']); nTrialsAborted=nTrialsAborted.nTrialsAborted;
    EoR_optimalities = load([mouseDataDir '/EoR_optimalities.mat']); EoR_optimalities=EoR_optimalities.EoR_optimalities;
    percentCompletedPR_allTrials = load([mouseDataDir '/percentCompletedPR_allTrials.mat']); percentCompletedPR_allTrials = percentCompletedPR_allTrials.percentCompletedPR_allTrials;
else
    NL_observed = load('~/phd/lever_task/optimality/WT/NL_observed.mat'); NL_observed=NL_observed.NL_observed;
    nTrialsAborted = load('~/phd/lever_task/optimality/WT/nTrialsAborted.mat'); nTrialsAborted=nTrialsAborted.nTrialsAborted;
    EoR_optimalities = load('~/phd/lever_task/optimality/WT/EoR_optimalities.mat'); EoR_optimalities=EoR_optimalities.EoR_optimalities;
    percentCompletedPR_allTrials = load('~/phd/lever_task/driftRL/abort_analysis/percentCompletedPR_allTrials.mat'); percentCompletedPR_allTrials = percentCompletedPR_allTrials.percentCompletedPR_allTrials;
end

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
            histogram(NL_observed{i},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
            hold on;
            if (p.Results.histograms_only)
                histogram(numLR{i},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
            else
                [F,xi,bw] = ksdensity(numLR{i},NUM_LR_VALS,'Bandwidth',NUM_LR_BANDWIDTH);
                plot(xi,F,'LineWidth',2)
            end
            xlabel('# of completed PR trials','FontSize',15,'FontWeight','bold')
            title([sessionTypes{i}])
        elseif (i > 4 & i <= 8)
            curInd = i - 4;
            histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
            hold on;
            if (p.Results.histograms_only)
                histogram(numAborted{curInd},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
            else
                [F,xi,bw] = ksdensity(numAborted{curInd},'Bandwidth',NUM_ABORTED_BANDWIDTH);
                plot(xi,F,'LineWidth',2)
            end
            xlabel('# of trials aborted','FontSize',15,'FontWeight','bold')
        elseif (i > 8 & i <= 12)
            curInd = i - 8;
            histogram(EoR_optimalities{curInd},'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2);
            hold on;
            if (p.Results.histograms_only)
                histogram(performanceEoR{curInd},'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2);
            else
                [F,xi,bw] = ksdensity(performanceEoR{curInd},PERFORMANCE_EOR_VALS,'Bandwidth',PERFORMANCE_EOR_BANDWIDTH);
                plot(xi,F,'LineWidth',2)
            end
            xlabel('EoR optimality','FontSize',15,'FontWeight','bold')
        else
            curInd = i - 12;
            histogram(percentCompletedPR_allTrials{curInd}(percentCompletedPR_allTrials{curInd} < 1),'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2)
            hold on;
            if (p.Results.histograms_only)
                histogram(percentCompletedPR{curInd}(percentCompletedPR{curInd} < 1),'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2)
            else
                [F,xi,bw] = ksdensity(percentCompletedPR{curInd}(percentCompletedPR{curInd} < 1),PERCENT_COMPLETED_PR_VALS,'Bandwidth',PERCENT_COMPLETED_PR_BANDWIDTH);
                plot(xi,F,'LineWidth',2)
            end
            xlabel('% of PR trial completed','FontSize',15,'FontWeight','bold')
        end
    end
    sgtitle(num2str(params))
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
                        histogram(NL_observed{j},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                        hold on;
                        if (p.Results.histograms_only)
                            histogram(numLR{j},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                        else
                            [F,xi,bw] = ksdensity(numLR{j},NUM_LR_VALS,'Bandwidth',NUM_LR_BANDWIDTH);
                            plot(xi,F,'LineWidth',2)
                        end
                        xlabel('# of completed PR trials','FontSize',15,'FontWeight','bold')
                        title([sessionTypes{j}])
                    elseif (j > 4 & j <= 8)
                        curInd = j - 4;
                        histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                        hold on;
                        if (p.Results.histograms_only)
                            histogram(numAborted{curInd},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                        else
                            [F,xi,bw] = ksdensity(numAborted{curInd},NUM_ABORTED_VALS,'Bandwidth',NUM_ABORTED_BANDWIDTH);
                            plot(xi,F,'LineWidth',2)
                        end
                        xlabel('# of trials aborted','FontSize',15,'FontWeight','bold')
                    else
                        curInd = j - 8;
                        histogram(EoR_optimalities{curInd},'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2);
                        hold on;
                        if (p.Results.histograms_only)
                            histogram(performanceEoR{curInd},'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2);
                        else
                            [F,xi,bw] = ksdensity(performanceEoR{curInd},PERFORMANCE_EOR_VALS,'Bandwidth',PERFORMANCE_EOR_BANDWIDTH);
                            plot(xi,F,'LineWidth',2)
                        end
                        xlabel('EoR optimality','FontSize',15,'FontWeight','bold')
                    end
                end
                curtitle = {num2str(params),['index = ' num2str(inds(i)) ' score = ' num2str(scores(inds(i)))]};
                sgtitle(curtitle)
                set(gcf,'Position',[10 10 1400 1200])
            elseif (strcmp(p.Results.fitMethod,'logprob_independent'))
                for j=1:16
                    subplot(4,4,j)
                    if (j <= 4)
                        histogram(NL_observed{j},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                        hold on;
                        if (p.Results.histograms_only)
                            histogram(numLR{j},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                        else
                            [F,xi,bw] = ksdensity(numLR{j},NUM_LR_VALS,'Bandwidth',NUM_LR_BANDWIDTH);
                            plot(xi,F,'LineWidth',2)
                        end
                        xlabel('# of completed PR trials','FontSize',15,'FontWeight','bold')
                        title([sessionTypes{j}])
                    elseif (j > 4 & j <= 8)
                        curInd = j - 4;
                        histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                        hold on;
                        if (p.Results.histograms_only)
                            histogram(numAborted{curInd},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                        else
                            [F,xi,bw] = ksdensity(numAborted{curInd},NUM_ABORTED_VALS,'Bandwidth',NUM_ABORTED_BANDWIDTH);
                            plot(xi,F,'LineWidth',2)
                        end
                        xlabel('# of trials aborted','FontSize',15,'FontWeight','bold')
                    elseif (j > 8 & j <= 12)
                        curInd = j - 8;
                        histogram(EoR_optimalities{curInd},'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2);
                        hold on;
                        if (p.Results.histograms_only)
                            histogram(performanceEoR{curInd},'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2);
                        else
                            [F,xi,bw] = ksdensity(performanceEoR{curInd},PERFORMANCE_EOR_VALS,'Bandwidth',PERFORMANCE_EOR_BANDWIDTH);
                            plot(xi,F,'LineWidth',2)
                        end
                        xlabel('EoR optimality','FontSize',15,'FontWeight','bold')
                    else
                        curInd = j - 12;
                        histogram(percentCompletedPR_allTrials{curInd}(percentCompletedPR_allTrials{curInd} < 1),'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2)
                        hold on;
                        if (p.Results.histograms_only)
                            histogram(percentCompletedPR{curInd}(percentCompletedPR{curInd} < 1),'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2)
                        else
                            [F,xi,bw] = ksdensity(percentCompletedPR{curInd}(percentCompletedPR{curInd} < 1),PERCENT_COMPLETED_PR_VALS,'Bandwidth',PERCENT_COMPLETED_PR_BANDWIDTH);
                            plot(xi,F,'LineWidth',2)
                        end
                        xlabel('% of PR trial completed','FontSize',15,'FontWeight','bold')
                    end
                end
                curtitle = {num2str(params),['index = ' num2str(inds(i)) ' score = ' num2str(scores(inds(i)))]};
                sgtitle(curtitle)
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
                    histogram(NL_observed{j},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                    hold on;
                    if (p.Results.histograms_only)
                        histogram(numLR{j},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                    else
                        [F,xi,bw] = ksdensity(numLR{j},NUM_LR_VALS,'Bandwidth',NUM_LR_BANDWIDTH);
                        plot(xi,F,'LineWidth',2)
                    end
                    xlabel('# of completed PR trials','FontSize',15,'FontWeight','bold')
                    title([sessionTypes{j}])
                elseif (j > 4 & j <= 8)
                    curInd = j - 4;
                    histogram(nTrialsAborted{curInd},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                    hold on;
                    if (p.Results.histograms_only)
                        histogram(numAborted{curInd},'normalization','pdf','binwidth',2,'DisplayStyle','stairs','LineWidth',2);
                    else
                        [F,xi,bw] = ksdensity(numAborted{curInd},NUM_ABORTED_VALS,'Bandwidth',NUM_ABORTED_BANDWIDTH);
                        plot(xi,F,'LineWidth',2)
                    end
                    xlabel('# of trials aborted','FontSize',15,'FontWeight','bold')
                elseif (j > 8 & j <= 12)
                    curInd = j - 8;
                    histogram(EoR_optimalities{curInd},'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2);
                    hold on;
                    if (p.Results.histograms_only)
                        histogram(performanceEoR{curInd},'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2);
                    else
                        [F,xi,bw] = ksdensity(performanceEoR{curInd},PERFORMANCE_EOR_VALS,'Bandwidth',PERFORMANCE_EOR_BANDWIDTH);
                        plot(xi,F,'LineWidth',2)
                    end
                    xlabel('EoR optimality','FontSize',15,'FontWeight','bold')
                else
                    curInd = j - 12;
                    histogram(percentCompletedPR_allTrials{curInd}(percentCompletedPR_allTrials{curInd} < 1),'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2)
                    hold on;
                    if (p.Results.histograms_only)
                        histogram(percentCompletedPR{curInd}(percentCompletedPR{curInd} < 1),'normalization','pdf','binwidth',.02,'DisplayStyle','stairs','LineWidth',2)
                    else
                        [F,xi,bw] = ksdensity(percentCompletedPR{curInd}(percentCompletedPR{curInd} < 1),PERCENT_COMPLETED_PR_VALS,'Bandwidth',PERCENT_COMPLETED_PR_BANDWIDTH);
                        plot(xi,F,'LineWidth',2)
                    end
                    xlabel('% of PR trial completed','FontSize',15,'FontWeight','bold')
                end
            end
            curtitle = {num2str(params),['index = ' num2str(inds(i)) ' score = ' num2str(scores(inds(i)))]};
            sgtitle(curtitle)
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

