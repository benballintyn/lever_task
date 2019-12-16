function [paramVals,paramNames,scores] = analzyeParamSlices(basedir,modelDir,subdir)
[paramVals,paramNames,scores] = getParamVals(basedir,modelDir,subdir);
info = load([basedir '/' modelDir '/' subdir '/' modelDir '.mat']); info=info.(modelDir);
npts = 2000;
[~,minInd] = min(scores);
for i=1:length(paramNames)
    for j=1:length(paramNames)
        if (i>=j)
            continue;
        else
            F = scatteredInterpolant(paramVals(:,i),paramVals(:,j),scores','natural');
            xmin = info.lower_bounds(i); xmax = info.upper_bounds(i);
            ymin = info.lower_bounds(j); ymax = info.upper_bounds(j);
            dx = (xmax-xmin)/npts;
            dy = (ymax-ymin)/npts;
            xrange=xmin:dx:xmax; yrange = ymin:dy:ymax;
            [Xq,Yq] = meshgrid(xrange,yrange);
            Vq = F(Xq,Yq);
            figure;
            imagesc([xmin xmax],[ymin ymax],Vq); colormap(jet); colorbar(); caxis([floor(min(scores)) ceil(max(scores))])
            set(gca,'ydir','normal')
            xlabel(paramNames{i}); ylabel(paramNames{j})
            hold on; plot(paramVals(:,i),paramVals(:,j),'.k','markersize',10)
            plot(paramVals(minInd,i),paramVals(minInd,j),'r*','markersize',20)
        end
    end
end
for i=1:size(paramVals,2)
    X(:,i) = paramVals(:,i);%/(info.upper_bounds(i) - info.lower_bounds(i));
end
%{
[minscore,minind] = min(scores);
mininds = find(scores == minscore);
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(X);
cmap=jet;
colorScale = ceil(max(scores))/size(cmap,1);
colorInds = ceil(scores/colorScale);
figure;
scatter3(paramVals*COEFF(:,1),paramVals*COEFF(:,2),paramVals*COEFF(:,3),[],cmap(colorInds,:),'filled')
colorScale = ceil((max(scores) - min(scores)))/size(cmap,1);
colorInds = ceil(((scores - min(scores)))/colorScale);
colorInds(mininds) = 1;
figure;
scatter3(paramVals*COEFF(:,1),paramVals*COEFF(:,2),paramVals*COEFF(:,3),[],cmap(colorInds,:),'filled')
hold on;
scatter3(paramVals(minind,:)*COEFF(:,1),paramVals(minind,:)*COEFF(:,2),paramVals(minind,:)*COEFF(:,3),200,'r*')
colormap(cmap); colorbar(); caxis([floor(min(scores)) ceil(max(scores))])
xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3')
figure;
plot(EXPLAINED)
title('Variance explained')
figure;
plot(COEFF)
title('PC coefficients')
legend(paramNames)
[~,inds] = sort(scores);
figure;
imagesc(X(inds,:))
set(gca,'xtick',1:size(paramVals,2),'xticklabel',paramNames)
colormap(jet); colorbar()
%}
meanParamVals = mean(paramVals,1);
stdParamVals = std(paramVals,[],1);
standardizedParamVals = (paramVals - meanParamVals)./stdParamVals;
standardizedScores = (scores - mean(scores))/std(scores);
[b,bint,r,rint,stats] = regress(standardizedScores',[standardizedParamVals ones(size(standardizedParamVals,1),1)]);

end

