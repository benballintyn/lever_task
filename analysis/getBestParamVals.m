function [C,meanScores] = getBestParamVals(basedir,modelDir,subdir,nClusters,N)
[paramVals,paramNames,scores,subscores] = getParamVals(basedir,modelDir,subdir);
[~,inds] = sort(scores);
Z = linkage(paramVals(inds(1:N),:));
dendrogram(Z)
[coeff,score,latent,tsquared,explained,mu] = pca(paramVals(inds(1:N),:));
figure;
scatter3(score(:,1),score(:,2),score(:,3),20,'filled')
[idx,C]=kmeans(paramVals(inds(1:N),:),nClusters);
for i=1:nClusters
    meanScores(i) = mean(scores(inds(idx==i)));
end
end

