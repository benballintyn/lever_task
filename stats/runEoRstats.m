function [isnormal,posthocPs] = runEoRstats(LOAD_PATH, SAVE_PATH)
eors = load([LOAD_PATH 'EoR_optimalities.mat']); eors=eors.EoR_optimalities;

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
save([SAVE_PATH 'p.mat'],'p','-mat')
save([SAVE_PATH 'H.mat'],'H','-mat')
save([SAVE_PATH 'p2.mat'],'p2','-mat')

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
            [p,h] = ranksum(eors{i},eors{j});
            posthocPs(i,j) = p*4;
        end
    end
end
posthocPs
save([SAVE_PATH 'posthocPs.mat'],'posthocPs','-mat')
end

