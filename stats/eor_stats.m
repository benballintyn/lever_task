% run stats on EoR optimalities

eors = load(['optimality/EoR_optimalities.mat']); eors=eors.EoR_optimalities;

for i=1:length(eors)
    [h,~] = adtest(eors{i});
    isnormal(i) = ~h;
end

count = 0;
for i=1:length(eors)
    for j=1:length(eors{i})
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
        Y(count) = eors{i}(j);
        LRs(count) = LR;
        Pss(count) = Ps;
    end
end

group = {LRs Pss};
[p,H,p2] = scheirer_ray_hare(Y,group)
save('stats/p.mat','p','-mat')
save('stats/H.mat','H','-mat')
save('stats/p2.mat','p2','-mat')

posthocPs = nan(4);
for i=1:4
    for j=1:4
        if (i==j)
            continue;
        else
            [p,h] = ranksum(eors{i},eors{j});
            posthocPs(i,j) = p*6;
        end
    end
end