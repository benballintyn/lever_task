function [] = findBestFitToMice(datadir,varargin)
isboolean = @(x) any(x == [0 1]);
p = inputParser;
addRequired(p,'datadir')
addParameter(p,'fitMethod','wasserstein',@ischar)
addParameter(p,'do_plot',false,isboolean)
addParameter(p,'redo',false,isboolean)
parse(p,datadir,varargin{:})

mice = {'AZ04','AZ08','HP01','HP02','HP03','HP04','MA01','NS07','NS09','NS10'};
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
NL_optimal = [10 22 28 58];
NL_observed = load('optimality/WT/NL_observed.mat'); NL_observed=NL_observed.NL_observed;
nL = load([datadir '/nL.mat']); nL=nL.nL;
numLR = load([datadir '/numLR.mat']); numLR=numLR.numLR;
performanceEoR = load([datadir '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;
performanceRoE = load([datadir '/performanceRoE.mat']); performanceRoE=performanceRoE.performanceRoE;
sweepParams = load([datadir '/sweepParams.mat']); sweepParams=sweepParams.sweepParams;
params = load([datadir '/params.mat']); params=params.params;

% Collect number of values per parameter into N
for i=1:length(params.var_names)
    N(i) = length(params.(params.var_names{i}));
end
N = [N 4]; % add another dimension with 4 entries for 4 session types

switch p.Results.fitMethod
    case 'wasserstein'
        for i=1:4
            [valcdf,valRange] = getCDF(NL_observed{i},1,1000,1); % create cdf of # of NL trials for the mice
            mouseCDF{i} = valcdf;
        end

        if (~exist([datadir '/bestMatchInds_wasserstein.mat'],'file') || p.Results.redo)
            % Compute wasserstein distances and find best parameter index values
            bestMatchInds_wasserstein = zeros(4,length(N)-1);
            bestEoRInds = zeros(4,length(N)-1);
            bestRoEInds = zeros(4,length(N)-1);
            wasserstein_distances = zeros(N);
            minwd = ones(1,4)*double(intmax);
            maxEoR = zeros(1,4);
            maxRoE = zeros(1,4);
            startTime=tic;
            nParamSets = prod(N);
            for i=1:nParamSets
                sessIndex = ceil(i/prod(N(1:end-1)));
                inds = cell(1,length(N));
                [inds{:}] = ind2sub(N,i);
                curInds = [inds{:}];
                [valcdf,valRange] = getCDF(numLR{i},1,1000,1);
                rlcdfs{i} = valcdf;
                wasserstein_distances(i) = wasserstein_1d(mouseCDF{sessIndex},valcdf);
                if (wasserstein_distances(i) < minwd(sessIndex))
                    minwd(sessIndex) = wasserstein_distances(i);
                    bestMatchInds_wasserstein(sessIndex,:) = curInds(1:end-1);
                end
                meanEoR = mean(performanceEoR{i});
                if (meanEoR > maxEoR(sessIndex))
                    maxEoR(sessIndex) = meanEoR;
                    bestEoRInds(sessIndex,:) = curInds(1:end-1);
                end
                meanRoE = mean(performanceRoE{i});
                if (meanRoE > maxRoE(sessIndex))
                    maxRoE(sessIndex) = meanRoE;
                    bestRoEInds(sessIndex,:) = curInds(1:end-1);
                end
                if (mod(i,1000) == 0)
                    elapsedTime = toc(startTime);
                    speed = elapsedTime/i;
                    timeRemaining_mins = speed*(nParamSets - i)/60;
                    disp(['Done with first ' num2str(i) ' parameter sets. ~' num2str(timeRemaining_mins) ' minutes left'])
                end
            end
            save([datadir '/wasserstein_distances.mat'],'wasserstein_distances','-mat')
            save([datadir '/bestMatchInds_wasserstein.mat'],'bestMatchInds_wasserstein','-mat')
            save([datadir '/bestEoRInds.mat'],'bestEoRInds','-mat')
            save([datadir '/bestRoEInds.mat'],'bestRoEInds','-mat')
        else
            disp('Wasserstein distances and best indices already found. Loading ...')
            wasserstein_distances = load([datadir '/wasserstein_distances.mat']); wasserstein_distances = wasserstein_distances.wasserstein_distances;
            bestMatchInds_wasserstein = load([datadir '/bestMatchInds_wasserstein.mat']); bestMatchInds_wasserstein = bestMatchInds_wasserstein.bestMatchInds_wasserstein;
            bestEoRInds = load([datadir '/bestEoRInds.mat']); bestEoRInds = bestEoRInds.bestEoRInds;
            bestRoEInds = load([datadir '/bestRoEInds.mat']); bestRoEInds = bestRoEInds.bestRoEInds;
        end
    case 'log_prob'
        if (~exist([datadir '/log_probs.mat']) || p.Results.redo)
            bestMatchInds_logProb = zeros(1,length(N)-1);
            bestEoRInds = zeros(4,length(N)-1);
            bestRoEInds = zeros(4,length(N)-1);
            log_probs = zeros(N);
            max_log_prob = -inf;
            maxEoR = zeros(1,4);
            maxRoE = zeros(1,4);
            startTime=tic;
            nParamSets = prod(N);
            for i=1:nParamSets
                sessIndex = ceil(i/prod(N(1:end-1)));
                inds = cell(1,length(N));
                [inds{:}] = ind2sub(N,i);
                curInds = [inds{:}];
                cur_log_prob = logprob2d(NL_observed{sessIndex},numLR{i});
                log_probs(i) = cur_log_prob;
                meanEoR = mean(performanceEoR{i});
                if (meanEoR > maxEoR(sessIndex))
                    maxEoR(sessIndex) = meanEoR;
                    bestEoRInds(sessIndex,:) = curInds(1:end-1);
                end
                meanRoE = mean(performanceRoE{i});
                if (meanRoE > maxRoE(sessIndex))
                    maxRoE(sessIndex) = meanRoE;
                    bestRoEInds(sessIndex,:) = curInds(1:end-1);
                end
                if (mod(i,1000) == 0)
                    elapsedTime = toc(startTime);
                    speed = elapsedTime/i;
                    timeRemaining_mins = speed*(nParamSets - i)/60;
                    disp(['Done with first ' num2str(i) ' parameter sets. ~' num2str(timeRemaining_mins) ' minutes left'])
                end
            end
            sumLogProbs = sum(log_probs,length(size(log_probs)));
            [~,bestInd] = max(sumLogProbs(:));
            inds = cell(1,length(N)-1);
            [inds{:}] = ind2sub(N(1:end-1),bestInd);
            bestMatchInds_logProb = [inds{:}];
            save([datadir '/log_probs.mat'],'log_probs','-mat')
            save([datadir '/bestMatchInds_logProb.mat'],'bestMatchInds_logProb','-mat')
            save([datadir '/bestEoRInds.mat'],'bestEoRInds','-mat')
            save([datadir '/bestRoEInds.mat'],'bestRoEInds','-mat')
        else
            disp('Log probabilities and best indices already found. Loading ...')
            log_probs = load([datadir '/log_probs.mat']); log_probs=log_probs.log_probs;
            bestMatchInds_logProb = load([datadir '/bestMatchInds_logProb.mat']); bestMatchInds_logProb = bestMatchInds_logProb.bestMatchInds_logProb;
            bestEoRInds = load([datadir '/bestEoRInds.mat']); bestEoRInds = bestEoRInds.bestEoRInds;
            bestRoEInds = load([datadir '/bestRoEInds.mat']); bestRoEInds = bestRoEInds.bestRoEInds;
        end
end



for i=1:length(params.var_names)
    if (strcmp(params.var_names{i},'ans_sigma'))
        adjustedVarNames{i} = 'ANS \sigma';
    else
        adjustedVarNames{i} = params.var_names{i};
    end
end
if (p.Results.do_plot)
    if (~exist([datadir '/figures'],'dir'))
        mkdir([datadir '/figures'])
    end
    switch p.Results.fitMethod
        case 'wasserstein'
            % Plot num PR histograms for real mice and best match parameters
            figure;
            for i=1:4
                subplot(2,2,i)
                linInd = sub2ind_cell(numLR,[bestMatchInds_wasserstein(i,:) i]);
                [h,p] = kstest2(NL_observed{i},numLR{linInd});
                histogram(NL_observed{i},'binwidth',2,'normalization','pdf')
                hold on;
                histogram(numLR{linInd},'binwidth',2,'normalization','pdf')
                yl = ylim;
                ylim([0 (yl(2) + yl(2)/5)])
                newyl = ylim;
                yval1 = yl(2) + (newyl(2) - yl(2))/2;
                yval2 = yl(2) + (newyl(2) - yl(2))/4;
                plot([(mean(NL_observed{i})-std(NL_observed{i})) (mean(numLR{linInd}) + std(numLR{linInd}))],[yval2 yval2],'Color','k')
                plot([NL_optimal(i) NL_optimal(i)],[0 newyl(2)],'Color','k')
                xval = mean([(mean(NL_observed{i})-std(NL_observed{i})) (mean(numLR{linInd}) + std(numLR{linInd}))]);
                if (p < .05 && p > .01)
                    text(xval,yval1,'*','fontsize',15)
                elseif (p < .01 && p > .001)
                    text(xval,yval1,'**','fontsize',15)
                elseif (p < .001)
                    text(xval,yval1,'***','fontsize',15)
                else
                    text(xval,yval1,'n.s.','fontsize',15)
                end
                xlim([0 100])
                % Create title str with each parameter value
                titleStr{1} = 'Best Match';
                for j=1:length(params.var_names)
                    titleStr{j+1} = [adjustedVarNames{j} ' = ' num2str(params.(params.var_names{j})(bestMatchInds_wasserstein(i,j)))];
                end
                title(titleStr,'fontsize',10,'fontweight','bold')
            end
            %suptitle('Best match to mice')
            set(gcf,'Position',[10 10 1400 1400])
            saveas(gcf,[datadir '/figures/num_PR_trials_best_match_wasserstein.fig'],'fig')
            saveas(gcf,[datadir '/figures/num_PR_trials_best_match_wasserstein.eps'],'eps')

            % Plot num PR histograms for real mice and best EoR parameters
            figure;
            for i=1:4
                subplot(2,2,i)
                histogram(NL_observed{i},'binwidth',2,'normalization','pdf');
                hold on;
                linInd = sub2ind_cell(numLR,[bestEoRInds(i,:) i]);
                histogram(numLR{linInd},'binwidth',2,'normalization','pdf')
                yl = ylim;
                plot([NL_optimal(i) NL_optimal(i)],[0 yl(2)],'Color','k')
                xlim([0 100])
                % Creat title str with each parameter value
                titleStr{1} = 'Best EoR';
                for j=1:length(params.var_names)
                    titleStr{j+1} = [adjustedVarNames{j} ' = ' num2str(params.(params.var_names{j})(bestEoRInds(i,j)))];
                end
                title(titleStr,'fontsize',10,'fontweight','bold')
            end
            %suptitle('Best EoR')
            set(gcf,'Position',[10 10 1400 1400])
            saveas(gcf,[datadir '/figures/num_PR_trials_bestEoRperformance.fig'],'fig')
            saveas(gcf,[datadir '/figures/num_PR_trials_bestEoRperformance.eps'],'eps')

            % plot num PR histograms for real mice and best RoE parameters
            figure;
            for i=1:4
                subplot(2,2,i)
                histogram(NL_observed{i},'binwidth',2,'normalization','pdf');
                hold on;
                linInd = sub2ind_cell(numLR,[bestRoEInds(i,:) i]);
                histogram(numLR{linInd},'binwidth',2,'normalization','pdf')
                yl = ylim;
                plot([NL_optimal(i) NL_optimal(i)],[0 yl(2)],'Color','k')
                xlim([0 100])
                % Creat title str with each parameter value
                titleStr{1} = 'Best RoE';
                for j=1:length(params.var_names)
                    titleStr{j+1} = [adjustedVarNames{j} ' = ' num2str(params.(params.var_names{j})(bestRoEInds(i,j)))];
                end
                title(titleStr,'fontsize',10,'fontweight','bold')
            end
            %suptitle('Best RoE')
            set(gcf,'Position',[10 10 1400 1400])
            saveas(gcf,[datadir '/figures/num_PR_trials_bestRoEperformance.fig'],'fig')
            saveas(gcf,[datadir '/figures/num_PR_trials_bestRoEperformance.eps'],'eps')

            % For each pair of parameters, plot the slice of wasserstein distances
            % for those 2 parameters keep all other parameters constant at their
            % best match values
            for i=1:length(params.var_names)
                for j=(i+1):length(params.var_names)
                    figure;
                    for k=1:4
                        S.subs = repmat({':'},1,ndims(wasserstein_distances));
                        unusedInds = setdiff(1:length(params.var_names),[i j]);
                        for l=1:length(unusedInds)
                            S.subs{unusedInds(l)} = bestMatchInds_wasserstein(k,unusedInds(l));
                        end
                        S.subs{end} = k;
                        S.type = '()';
                        wd_slice = squeeze(subsref(wasserstein_distances,S));
                        subplot(2,2,k)
                        surf(params.(params.var_names{i}),params.(params.var_names{j}),wd_slice')
                        xlabel(adjustedVarNames{i})
                        ylabel(adjustedVarNames{j})
                    end
                    suptitle('Best match to mice')
                    set(gcf,'Position',[10 10 1400 1400])
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestMatch_wasserstein_slice.fig'],'fig')
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestMatch_wasserstein_slice.eps'],'eps')
                end
            end

            % For each pair of parameters, plot the slice of wasserstein distances
            % for those 2 parameters keep all other parameters constant at their
            % best EoR values
            for i=1:length(params.var_names)
                for j=(i+1):length(params.var_names)
                    figure;
                    for k=1:4
                        S.subs = repmat({':'},1,ndims(wasserstein_distances));
                        unusedInds = setdiff(1:length(params.var_names),[i j]);
                        for l=1:length(unusedInds)
                            S.subs{unusedInds(l)} = bestEoRInds(k,unusedInds(l));
                        end
                        S.subs{end} = k;
                        S.type = '()';
                        wd_slice = squeeze(subsref(wasserstein_distances,S));
                        subplot(2,2,k)
                        surf(params.(params.var_names{i}),params.(params.var_names{j}),wd_slice')
                        xlabel(adjustedVarNames{i})
                        ylabel(adjustedVarNames{j})
                    end
                    suptitle('Best EoR')
                    set(gcf,'Position',[10 10 1400 1400])
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestEoR_wasserstein_slice.fig'],'fig')
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestEoR_wasserstein_slice.eps'],'eps')
                end
            end

            % For each pair of parameters, plot the slice of wasserstein distances
            % for those 2 parameters keep all other parameters constant at their
            % best RoE values
            for i=1:length(params.var_names)
                for j=(i+1):length(params.var_names)
                    S.subs = repmat({':'},1,ndims(wasserstein_distances));
                    unusedInds = setdiff(1:length(params.var_names),[i j]);
                    figure;
                    for k=1:4
                        S.subs = repmat({':'},1,ndims(wasserstein_distances));
                        unusedInds = setdiff(1:length(params.var_names),[i j]);
                        for l=1:length(unusedInds)
                            S.subs{unusedInds(l)} = bestRoEInds(k,unusedInds(l));
                        end
                        S.subs{end} = k;
                        S.type = '()';
                        wd_slice = squeeze(subsref(wasserstein_distances,S));
                        subplot(2,2,k)
                        surf(params.(params.var_names{i}),params.(params.var_names{j}),wd_slice')
                        xlabel(adjustedVarNames{i})
                        ylabel(adjustedVarNames{j})
                    end
                    suptitle('Best RoE')
                    set(gcf,'Position',[10 10 1400 1400])
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestRoE_wasserstein_slice.fig'],'fig')
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestRoE_wasserstein_slice.eps'],'eps')
                end
            end
        case 'log_prob'
            % Plot num PR histograms for real mice and best match parameters
            figure;
            for i=1:4
                subplot(2,2,i)
                linInd = sub2ind_cell(numLR,[bestMatchInds_logProb i]);
                [h,p] = kstest2(NL_observed{i},numLR{linInd});
                histogram(NL_observed{i},'binwidth',2,'normalization','pdf')
                hold on;
                histogram(numLR{linInd},'binwidth',2,'normalization','pdf')
                yl = ylim;
                ylim([0 (yl(2) + yl(2)/5)])
                newyl = ylim;
                yval1 = yl(2) + (newyl(2) - yl(2))/2;
                yval2 = yl(2) + (newyl(2) - yl(2))/4;
                plot([(mean(NL_observed{i})-std(NL_observed{i})) (mean(numLR{linInd}) + std(numLR{linInd}))],[yval2 yval2],'Color','k')
                plot([NL_optimal(i) NL_optimal(i)],[0 newyl(2)],'Color','k')
                xval = mean([(mean(NL_observed{i})-std(NL_observed{i})) (mean(numLR{linInd}) + std(numLR{linInd}))]);
                if (p < .05 && p > .01)
                    text(xval,yval1,'*','fontsize',15)
                elseif (p < .01 && p > .001)
                    text(xval,yval1,'**','fontsize',15)
                elseif (p < .001)
                    text(xval,yval1,'***','fontsize',15)
                else
                    text(xval,yval1,'n.s.','fontsize',15)
                end
                xlim([0 100])
                % Create title str with each parameter value
                titleStr{1} = 'Best Match';
                for j=1:length(params.var_names)
                    titleStr{j+1} = [adjustedVarNames{j} ' = ' num2str(params.(params.var_names{j})(bestMatchInds_logProb(j)))];
                end
                title(titleStr,'fontsize',10,'fontweight','bold')
            end
            %suptitle('Best match to mice')
            set(gcf,'Position',[10 10 1400 1400])
            saveas(gcf,[datadir '/figures/num_PR_trials_best_match_logProb.fig'],'fig')
            saveas(gcf,[datadir '/figures/num_PR_trials_best_match_logProb.eps'],'eps')
            
            % Plot num PR histograms for real mice and best EoR parameters
            figure;
            for i=1:4
                subplot(2,2,i)
                histogram(NL_observed{i},'binwidth',2,'normalization','pdf');
                hold on;
                linInd = sub2ind_cell(numLR,[bestEoRInds(i,:) i]);
                histogram(numLR{linInd},'binwidth',2,'normalization','pdf')
                yl = ylim;
                plot([NL_optimal(i) NL_optimal(i)],[0 yl(2)],'Color','k')
                xlim([0 100])
                % Creat title str with each parameter value
                titleStr{1} = 'Best EoR';
                for j=1:length(params.var_names)
                    titleStr{j+1} = [adjustedVarNames{j} ' = ' num2str(params.(params.var_names{j})(bestEoRInds(i,j)))];
                end
                title(titleStr,'fontsize',10,'fontweight','bold')
            end
            %suptitle('Best EoR')
            set(gcf,'Position',[10 10 1400 1400])
            saveas(gcf,[datadir '/figures/num_PR_trials_bestEoRperformance.fig'],'fig')
            saveas(gcf,[datadir '/figures/num_PR_trials_bestEoRperformance.eps'],'eps')

            % plot num PR histograms for real mice and best RoE parameters
            figure;
            for i=1:4
                subplot(2,2,i)
                histogram(NL_observed{i},'binwidth',2,'normalization','pdf');
                hold on;
                linInd = sub2ind_cell(numLR,[bestRoEInds(i,:) i]);
                histogram(numLR{linInd},'binwidth',2,'normalization','pdf')
                yl = ylim;
                plot([NL_optimal(i) NL_optimal(i)],[0 yl(2)],'Color','k')
                xlim([0 100])
                % Creat title str with each parameter value
                titleStr{1} = 'Best RoE';
                for j=1:length(params.var_names)
                    titleStr{j+1} = [adjustedVarNames{j} ' = ' num2str(params.(params.var_names{j})(bestRoEInds(i,j)))];
                end
                title(titleStr,'fontsize',10,'fontweight','bold')
            end
            %suptitle('Best RoE')
            set(gcf,'Position',[10 10 1400 1400])
            saveas(gcf,[datadir '/figures/num_PR_trials_bestRoEperformance.fig'],'fig')
            saveas(gcf,[datadir '/figures/num_PR_trials_bestRoEperformance.eps'],'eps')
            
            % For each pair of parameters, plot the slice of log probabilities
            % for those 2 parameters keep all other parameters constant at their
            % best match values
            for i=1:length(params.var_names)
                for j=(i+1):length(params.var_names)
                    figure;
                    for k=1:4
                        S.subs = repmat({':'},1,ndims(log_probs));
                        unusedInds = setdiff(1:length(params.var_names),[i j]);
                        for l=1:length(unusedInds)
                            S.subs{unusedInds(l)} = bestMatchInds_logProb(unusedInds(l));
                        end
                        S.subs{end} = k;
                        S.type = '()';
                        lp_slice = squeeze(subsref(log_probs,S));
                        subplot(2,2,k)
                        surf(params.(params.var_names{i}),params.(params.var_names{j}),lp_slice')
                        xlabel(adjustedVarNames{i})
                        ylabel(adjustedVarNames{j})
                    end
                    suptitle('Best match to mice. Log prob')
                    set(gcf,'Position',[10 10 1400 1400])
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestMatch_logProb_slice.fig'],'fig')
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestMatch_logProb_slice.eps'],'eps')
                end
            end
            
            % For each pair of parameters, plot the slice of log probabilities
            % for those 2 parameters keep all other parameters constant at their
            % best EoR values
            for i=1:length(params.var_names)
                for j=(i+1):length(params.var_names)
                    figure;
                    for k=1:4
                        S.subs = repmat({':'},1,ndims(log_probs));
                        unusedInds = setdiff(1:length(params.var_names),[i j]);
                        for l=1:length(unusedInds)
                            S.subs{unusedInds(l)} = bestEoRInds(k,unusedInds(l));
                        end
                        S.subs{end} = k;
                        S.type = '()';
                        lp_slice = squeeze(subsref(log_probs,S));
                        subplot(2,2,k)
                        surf(params.(params.var_names{i}),params.(params.var_names{j}),lp_slice')
                        xlabel(adjustedVarNames{i})
                        ylabel(adjustedVarNames{j})
                    end
                    suptitle('Best EoR')
                    set(gcf,'Position',[10 10 1400 1400])
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestEoR_logProb_slice.fig'],'fig')
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestEoR_logProb_slice.eps'],'eps')
                end
            end
            
            % For each pair of parameters, plot the slice of log probabilities
            % for those 2 parameters keep all other parameters constant at their
            % best RoE values
            for i=1:length(params.var_names)
                for j=(i+1):length(params.var_names)
                    S.subs = repmat({':'},1,ndims(log_probs));
                    unusedInds = setdiff(1:length(params.var_names),[i j]);
                    figure;
                    for k=1:4
                        S.subs = repmat({':'},1,ndims(log_probs));
                        unusedInds = setdiff(1:length(params.var_names),[i j]);
                        for l=1:length(unusedInds)
                            S.subs{unusedInds(l)} = bestRoEInds(k,unusedInds(l));
                        end
                        S.subs{end} = k;
                        S.type = '()';
                        lp_slice = squeeze(subsref(log_probs,S));
                        subplot(2,2,k)
                        surf(params.(params.var_names{i}),params.(params.var_names{j}),lp_slice')
                        xlabel(adjustedVarNames{i})
                        ylabel(adjustedVarNames{j})
                    end
                    suptitle('Best RoE')
                    set(gcf,'Position',[10 10 1400 1400])
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestRoE_logProb_slice.fig'],'fig')
                    saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestRoE_logProb_slice.eps'],'eps')
                end
            end
    end
    
    
end
end

