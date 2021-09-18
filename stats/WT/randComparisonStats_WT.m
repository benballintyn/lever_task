% randComparisonStats
EoRrandomComparison = load(['~/phd/lever_task/optimality/WT/EoRrandomComparison.mat']); EoRrandomComparison=EoRrandomComparison.EoRrandomComparison;

for i=1:length(EoRrandomComparison)
    [h,~] = adtest(EoRrandomComparison{i});
    isnormal(i) = ~h;
end

count = 0;
for i=1:length(EoRrandomComparison)
    for j=1:length(EoRrandomComparison{i})
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
        Y(count) = EoRrandomComparison{i}(j);
        LRs(count) = LR;
        Pss(count) = Ps;
    end
end

group = {LRs Pss};
[p,H,p2] = scheirer_ray_hare(Y,group)
save('stats/p_randComparison.mat','p','-mat')
save('stats/H_randComparison.mat','H','-mat')
save('stats/p2_randComparison.mat','p2','-mat')

posthocPs = nan(4);
for i=1:4
    for j=1:4
        if (i==j)
            continue;
        elseif ((i==1 && j == 4) || (i == 4 && j == 1))
            posthocPs(i,j) = nan;
        elseif ((i==2 && j == 3) || (i == 3 && j == 2))
            posthocPs(i,j) = nan;
        else
            [p,h] = ranksum(EoRrandomComparison{i},EoRrandomComparison{j});
            posthocPs(i,j) = p*4;
        end
    end
end
posthocPs
save('stats/posthocPs_randComparison.mat','posthocPs','-mat')