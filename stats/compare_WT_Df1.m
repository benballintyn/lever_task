% compare_WT_Df1
clear all;
savedir = '~/phd/lever_task/stats/WT_Df1_comparison/';
if (~exist(savedir,'dir'))
    mkdir(savedir)
end
eorsWT = load(['optimality/WT/EoR_optimalities.mat']); eorsWT=eorsWT.EoR_optimalities;
eorsDf1 = load(['optimality/Df1/EoR_optimalities.mat']); eorsDf1=eorsDf1.EoR_optimalities;

EoRrandomComparisonWT = load(['~/phd/lever_task/optimality/WT/EoRrandomComparison.mat']); EoRrandomComparisonWT=EoRrandomComparisonWT.EoRrandomComparison;
EoRrandomComparisonDf1 = load(['~/phd/lever_task/optimality/Df1/EoRrandomComparison.mat']); EoRrandomComparisonDf1=EoRrandomComparisonDf1.EoRrandomComparison;

count = 0;
for i=1:length(eorsWT)
    for j=1:length(eorsWT{i})
        count=count+1;
        switch i
            case 1
                LR = 6;
                Ps = 6;
            case 2
                LR = 6;
                Ps = 12;
            case 3
                LR = 15;
                Ps = 6;
            case 4
                LR = 15;
                Ps = 12;
        end
        Y(count) = eorsWT{i}(j);
        Y2(count) = EoRrandomComparisonWT{i}(j);
        LRs(count) = LR;
        Pss(count) = Ps;
        mouseType{count} = 'WT';
    end
end

for i=1:length(eorsDf1)
    for j=1:length(eorsDf1{i})
        count=count+1;
        switch i
            case 1
                LR = 6;
                Ps = 6;
            case 2
                LR = 6;
                Ps = 12;
            case 3
                LR = 15;
                Ps = 6;
            case 4
                LR = 15;
                Ps = 12;
        end
        Y(count) = eorsDf1{i}(j);
        Y2(count) = EoRrandomComparisonDf1{i}(j);
        LRs(count) = LR;
        Pss(count) = Ps;
        mouseType{count} = 'Df1';
    end
end

group = {LRs Pss mouseType};
[p,H,p2] = scheirer_ray_hare(Y,group)
save([savedir 'p.mat'],'p','-mat')
save([savedir 'H.mat'],'H','-mat')
save([savedir 'p2.mat'],'p2','-mat')

[p,H,p2] = scheirer_ray_hare(Y2,group)
save([savedir 'p_randomComparison.mat'],'p','-mat')
save([savedir 'H_randomComparison.mat'],'H','-mat')
save([savedir 'p2_randomComparison.mat'],'p2','-mat')

save([savedir 'Y.mat'],'Y','-mat')
save([savedir 'Y2.mat'],'Y2','-mat')
save([savedir 'group.mat'],'group','-mat')

