function [score,subscores] = fake_data_objective_score(fakeDataDir,modelDir)
counts3D = load([fakeDataDir '/counts3D.mat']); counts3D=counts3D.counts3D;

numLR = load([modelDir '/numLR.mat']); numLR=numLR.numLR;
numAborted = load([modelDir '/numAborted.mat']); numAborted=numAborted.numAborted;
performanceEoR = load([modelDir '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;

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
    if (any(size(counts3D{i}) ~= size(modelData3D{i})))
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

    nonZeroInds = find(counts3D{i} ~= 0);
    nonZeroCounts = counts3D{i}(nonZeroInds);
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

