function [score,subscores] = rl_objective_score(datadir,scoreType,varargin)
p = inputParser;
addRequired(p,'datadir',@ischar)
addRequired(p,'scoreType',@ischar)
addParameter(p,'Only120Trials',false,@islogical)
parse(p,datadir,scoreType,varargin{:})

switch scoreType
    case 'logprob_independent'
        if (p.Results.Only120Trials)
            NL_observed = load('~/phd/lever_task/analysis_120_trials/NL_observed.mat'); NL_observed=NL_observed.NL_observed;
            nTrialsAborted = load('optimality/WT/nTrialsAborted.mat'); nTrialsAborted=nTrialsAborted.nTrialsAborted;
            EoR_optimalities = load('optimality/WT/EoR_optimalities.mat'); EoR_optimalities=EoR_optimalities.EoR_optimalities;
            %percentCompletedPR = load('driftRL/abort_analysis/percentCompletedPR.mat'); percentCompletedPR=percentCompletedPR.percentCompletedPR;
            percentCompletedPR_allTrials = load('driftRL/abort_analysis/percentCompletedPR_allTrials.mat'); percentCompletedPR_allTrials = percentCompletedPR_allTrials.percentCompletedPR_allTrials;
        else
            NL_observed = load('optimality/WT/NL_observed.mat'); NL_observed=NL_observed.NL_observed;
            nTrialsAborted = load('optimality/WT/nTrialsAborted.mat'); nTrialsAborted=nTrialsAborted.nTrialsAborted;
            EoR_optimalities = load('optimality/WT/EoR_optimalities.mat'); EoR_optimalities=EoR_optimalities.EoR_optimalities;
            %percentCompletedPR = load('driftRL/abort_analysis/percentCompletedPR.mat'); percentCompletedPR=percentCompletedPR.percentCompletedPR;
            percentCompletedPR_allTrials = load('driftRL/abort_analysis/percentCompletedPR_allTrials.mat'); percentCompletedPR_allTrials = percentCompletedPR_allTrials.percentCompletedPR_allTrials;
        end

        numLR = load([datadir '/numLR.mat']); numLR=numLR.numLR;
        numAborted = load([datadir '/numAborted.mat']); numAborted=numAborted.numAborted;
        performanceEoR = load([datadir '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;
        agentPercentCompletedPR = load([datadir '/percentCompletedPR.mat']); agentPercentCompletedPR=agentPercentCompletedPR.percentCompletedPR;

        logprobs_numLR = zeros(1,4);
        logprobs_aborted = zeros(1,4);
        logprobs_eor = zeros(1,4);
        logprobs_percentCompletedPR = zeros(1,4);
        for i=1:4
            [curlogprob_numLR] = logprob2d(NL_observed{i},numLR{i},'bandwidth',1);
            logprobs_numLR(i) = curlogprob_numLR;

            [curlogprob_aborted] = logprob2d(nTrialsAborted{i},numAborted{i},'bandwidth',2);
            logprobs_aborted(i) = curlogprob_aborted;

            [curlogprob_eor] = logprob2d(EoR_optimalities{i},performanceEoR{i},'bandwidth',.01);
            logprobs_eor(i) = curlogprob_eor;

            [curlogprob_percentCompletedPR] = logprob2d(percentCompletedPR_allTrials{i},agentPercentCompletedPR{i},'bandwidth',.01);
            logprobs_percentCompletedPR(i) = curlogprob_percentCompletedPR;
        end
        subscores(1,:) = logprobs_numLR;
        subscores(2,:) = logprobs_aborted;
        subscores(3,:) = logprobs_eor;
        subscores(4,:) = logprobs_percentCompletedPR;

        score = -sum(subscores(:));
        
    case 'logprob_joint'
        % Load mouse and model data
        if (p.Results.Only120Trials)
            mouseData3D_counts = load('/home/ben/phd/lever_task/analysis_120_trials/mouseData3D_counts_120Trials.mat');
            mouseData3D_counts = mouseData3D_counts.mouseData3D_counts_120Trials;
        else
            mouseData3D_counts = load('/home/ben/phd/lever_task/driftRL/mouseData3D_counts.mat'); 
            mouseData3D_counts=mouseData3D_counts.mouseData3D_counts;
        end
        
        numLR = load([datadir '/numLR.mat']); numLR=numLR.numLR;
        numAborted = load([datadir '/numAborted.mat']); numAborted=numAborted.numAborted;
        performanceEoR = load([datadir '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;
        
        if (p.Results.Only120Trials)
            nPR_counts = 15; nPR_binSize = 5;
            nAborted_counts = 14; nAborted_binSize = 5;
            EoR_counts = 20; EoR_binSize = .025;
            % Create bins
            % N_PR
            bins_lb{1} = 0:nPR_binSize:(nPR_binSize*(nPR_counts - 1));
            bins_ub{1} = nPR_binSize:nPR_binSize:(nPR_binSize*nPR_counts);
            
            % N_aborted
            bins_lb{2} = 0:nAborted_binSize:(nAborted_binSize*(nAborted_counts - 1));
            bins_ub{2} = nAborted_binSize:nAborted_binSize:(nAborted_binSize*nAborted_counts);
            
            % EoR optimality
            bins_lb{3} = .5:EoR_binSize:(.5 + (EoR_binSize*(EoR_counts - 1)));
            bins_ub{3} = (.5 + EoR_binSize):EoR_binSize:(.5 + EoR_binSize*EoR_counts);
            
            modelData3D = cell(1,4);
            logprobs = zeros(1,4);
            for i=1:4
                modelData3D{i} = zeros(length(bins_lb{1}),length(bins_lb{2}),length(bins_lb{3}));
                if (any(size(mouseData3D_counts{i}) ~= size(modelData3D{i})))
                    error('Model distribution is of a different size than the mouse distribution')
                end
                for j=1:length(numLR{i})
                    binInd = intersect(find(bins_lb{1} <= numLR{i}(j)),find(bins_ub{1} > numLR{i}(j)));
                    if (~isempty(binInd))
                        model_bin1{i}(j) = binInd;
                    else
                        model_bin1{i}(j) = length(bins_lb{1});
                    end
                    binInd = intersect(find(bins_lb{2} <= numAborted{i}(j)),find(bins_ub{2} > numAborted{i}(j)));
                    if (~isempty(binInd))
                        model_bin2{i}(j) = binInd;
                    else
                        model_bin2{i}(j) = length(bins_lb{2});
                    end
                    binInd = intersect(find(bins_lb{3} <= performanceEoR{i}(j)),find(bins_ub{3} > performanceEoR{i}(j)));
                    if (~isempty(binInd))
                        model_bin3{i}(j) = binInd;
                    else
                        model_bin3{i}(j) = length(bins_lb{3});
                    end
                    modelData3D{i}(model_bin1{i}(j),model_bin2{i}(j),model_bin3{i}(j)) = modelData3D{i}(model_bin1{i}(j),model_bin2{i}(j),model_bin3{i}(j)) + 1;
                end
                modelData3D{i} = modelData3D{i}/length(numLR{i});
                modelData3D{i}(modelData3D{i} == 0) = 1/realmax;
                modelData3D{i} = modelData3D{i}./sum(modelData3D{i}(:));

                nonZeroInds = find(mouseData3D_counts{i} ~= 0);
                nonZeroCounts = mouseData3D_counts{i}(nonZeroInds);
                probs = [];
                for j=1:length(nonZeroCounts)
                    for k=1:nonZeroCounts(j)
                        probs = [probs modelData3D{i}(nonZeroInds(j))];
                    end
                end
                logprobs(i) = sum(log(probs));
            end
            subscores = logprobs;
            score = -sum(logprobs);
        else
            % Create bins
            % N_PR
            bins_lb{1} = 0:5:95;
            bins_ub{1} = 5:5:100;

            % N_aborted
            bins_lb{2} = 0:20:180;
            bins_ub{2} = 20:20:200;

            % EoR optimality
            bins_lb{3} = 0:.05:.95;
            bins_ub{3} = .05:.05:1;

            modelData3D = cell(1,4);
            logprobs = zeros(1,4);
            for i=1:4
                modelData3D{i} = zeros(length(bins_lb{1}),length(bins_lb{2}),length(bins_lb{3}));
                for j=1:length(numLR{i})
                    binInd = intersect(find(bins_lb{1} <= numLR{i}(j)),find(bins_ub{1} > numLR{i}(j)));
                    if (~isempty(binInd))
                        model_bin1{i}(j) = binInd;
                    else
                        model_bin1{i}(j) = length(bins_lb{1});
                    end
                    binInd = intersect(find(bins_lb{2} <= numAborted{i}(j)),find(bins_ub{2} > numAborted{i}(j)));
                    if (~isempty(binInd))
                        model_bin2{i}(j) = binInd;
                    else
                        model_bin2{i}(j) = length(bins_lb{2});
                    end
                    binInd = intersect(find(bins_lb{3} <= performanceEoR{i}(j)),find(bins_ub{3} > performanceEoR{i}(j)));
                    if (~isempty(binInd))
                        model_bin3{i}(j) = binInd;
                    else
                        model_bin3{i}(j) = length(bins_lb{3});
                    end
                    modelData3D{i}(model_bin1{i}(j),model_bin2{i}(j),model_bin3{i}(j)) = modelData3D{i}(model_bin1{i}(j),model_bin2{i}(j),model_bin3{i}(j)) + 1;
                end
                modelData3D{i} = modelData3D{i}/length(numLR{i});
                modelData3D{i}(modelData3D{i} == 0) = 1/realmax;
                modelData3D{i} = modelData3D{i}./sum(modelData3D{i}(:));

                nonZeroInds = find(mouseData3D_counts{i} ~= 0);
                nonZeroCounts = mouseData3D_counts{i}(nonZeroInds);
                probs = [];
                for j=1:length(nonZeroCounts)
                    for k=1:nonZeroCounts(j)
                        probs = [probs modelData3D{i}(nonZeroInds(j))];
                    end
                end
                logprobs(i) = sum(log(probs));
            end
            subscores = logprobs;
            score = -sum(logprobs);
        end
end

