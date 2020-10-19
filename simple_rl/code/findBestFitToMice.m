function [] = findBestFitToMice(datadir,varargin)
isboolean = @(x) any(x == [0 1]);
p = inputParser;
addRequired(p,'datadir')
addParameter(p,'do_plot',false,isboolean)
addParameter(p,'redo',false,isboolean)
parse(p,datadir,varargin{:})

mice = {'AZ04','AZ08','HP01','HP02','HP03','HP04','MA01','NS07','NS09','NS10'};
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
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

for i=1:4
    [valcdf,valRange] = getCDF(NL_observed{i},1,1000,1); % create cdf of # of NL trials for the mice
    mouseCDF{i} = valcdf;
end

if (~exist([datadir '/wasserstein_distances.mat'],'file') || p.Results.redo)
    % Compute wasserstein distances and find best parameter index values
    bestMatchInds = zeros(4,length(N)-1);
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
            bestMatchInds(sessIndex,:) = curInds(1:end-1);
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
    save([datadir '/bestMatchInds.mat'],'bestMatchInds','-mat')
    save([datadir '/bestEoRInds.mat'],'bestEoRInds','-mat')
    save([datadir '/bestRoEInds.mat'],'bestRoEInds','-mat')
else
    disp('Wasserstein distances and best indices already found. Loading ...')
    wasserstein_distances = load([datadir '/wasserstein_distances.mat']); wasserstein_distances = wasserstein_distances.wasserstein_distances;
    bestMatchInds = load([datadir '/bestMatchInds.mat']); bestMatchInds = bestMatchInds.bestMatchInds;
    bestEoRInds = load([datadir '/bestEoRInds.mat']); bestEoRInds = bestEoRInds.bestEoRInds;
    bestRoEInds = load([datadir '/bestRoEInds.mat']); bestRoEInds = bestRoEInds.bestRoEInds;
end

if (p.Results.do_plot)
    if (~exist([datadir '/figures'],'dir'))
        mkdir([datadir '/figures'])
    end
    
    % Plot num PR histograms for real mice and best match parameters
    figure;
    for i=1:4
        subplot(2,2,i)
        histogram(NL_observed{i},'binwidth',2,'normalization','pdf')
        hold on;
        linInd = sub2ind_cell(numLR,[bestMatchInds(i,:) i]);
        histogram(numLR{linInd},'binwidth',2,'normalization','pdf')
        % Create title str with each parameter value
        titleStr{1} = 'Best Match';
        for j=1:length(params.var_names)
            titleStr{j+1} = [params.var_names{j} ' = ' num2str(params.(params.var_names{j})(bestMatchInds(i,j)))];
        end
        title(titleStr,'fontsize',10,'fontweight','bold')
    end
    %suptitle('Best match to mice')
    set(gcf,'Position',[10 10 1400 1400])
    saveas(gcf,[datadir '/figures/num_PR_trials_best_match.fig'],'fig')
    saveas(gcf,[datadir '/figures/num_PR_trials_best_match.eps'],'eps')
    
    % Plot num PR histograms for real mice and best EoR parameters
    figure;
    for i=1:4
        subplot(2,2,i)
        histogram(NL_observed{i},'binwidth',2,'normalization','pdf');
        hold on;
        linInd = sub2ind_cell(numLR,[bestEoRInds(i,:) i]);
        histogram(numLR{linInd},'binwidth',2,'normalization','pdf')
        % Creat title str with each parameter value
        titleStr{1} = 'Best EoR';
        for j=1:length(params.var_names)
            titleStr{j+1} = [params.var_names{j} ' = ' num2str(params.(params.var_names{j})(bestEoRInds(i,j)))];
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
        % Creat title str with each parameter value
        titleStr{1} = 'Best RoE';
        for j=1:length(params.var_names)
            titleStr{j+1} = [params.var_names{j} ' = ' num2str(params.(params.var_names{j})(bestRoEInds(i,j)))];
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
                    S.subs{unusedInds(l)} = bestMatchInds(k,unusedInds(l));
                end
                S.subs{end} = k;
                S.type = '()';
                wd_slice = squeeze(subsref(wasserstein_distances,S));
                subplot(2,2,k)
                surf(params.(params.var_names{i}),params.(params.var_names{j}),wd_slice')
                xlabel(params.var_names{i})
                ylabel(params.var_names{j})
            end
            suptitle('Best match to mice')
            set(gcf,'Position',[10 10 1400 1400])
            saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestMatch_slice.fig'],'fig')
            saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestMatch_slice.eps'],'eps')
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
                xlabel(params.var_names{i})
                ylabel(params.var_names{j})
            end
            suptitle('Best EoR')
            set(gcf,'Position',[10 10 1400 1400])
            saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestEoR_slice.fig'],'fig')
            saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestEoR_slice.eps'],'eps')
        end
    end
    
    % For each pair of parameters, plot the slice of wasserstein distances
    % for those 2 parameters keep all other parameters constant at their
    % best EoR values
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
                xlabel(params.var_names{i})
                ylabel(params.var_names{j})
            end
            suptitle('Best RoE')
            set(gcf,'Position',[10 10 1400 1400])
            saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestRoE_slice.fig'],'fig')
            saveas(gcf,[datadir '/figures/' params.var_names{i} '_' params.var_names{j} '_bestRoE_slice.eps'],'eps')
        end
    end
    
    
end
end

