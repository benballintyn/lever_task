% pLpS_crosses
pL = load('optimality/pL.mat'); pL=pL.pL;
pS = load('optimality/pS.mat'); pS=pS.pS;

for i=1:4
    smoothPL{i} = movmean(pL{i},3);
    smoothPS{i} = movmean(pS{i},3);
end

for i=1:4
    z = find(smoothPS{i}(1:200) < smoothPL{i}(1:200));
    cx(i) = z(end);
end
plot(cx)

figure;
for i=1:4
    subplot(2,2,i)
    plot(smoothPS{i}(1:500)); hold on; plot(smoothPL{i}(1:500))
    scatter(cx(i),.5,40,'ko','filled')
    set(gca,'xtick',[1 100 200 300 400 500],'xticklabels',{'1','','','','','500'});
    set(gca,'ytick',[0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1],'yticklabels',{'0','','','','','','','','','','1'})
    if (i==1)
        legend({'P(SR)', 'P(LR)'})
    end
    xlabel('Trial #')
end