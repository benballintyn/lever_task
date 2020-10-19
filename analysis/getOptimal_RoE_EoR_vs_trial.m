function [roe,eor] = getOptimal_RoE_EoR_vs_trial(Ntrials)
f_NL = @(NL) .5*(sqrt(8*NL+9) - 3);
RoE_rstar = @(N,NL,SR,LR,Ps) (N - NL)*(SR/Ps) + LR*(.5*(sqrt(8*NL+9)-3));
RoE_NL_optimal = [71 287 449 1799];
EoR_NL_optimal = [10 22 28 58];
eg = double(vpa(eulergamma));
EoR_star = @(N,NL,SR,LR,Ps) (1/N)*(LR*(psi(0,NL+2) + eg - 1) + (N - NL)*(SR/Ps));

for i=1:4
    if (i==1)
        SR = 3;
        LR = 6;
        Ps = 6;
    elseif (i==2)
        SR = 3;
        LR = 6;
        Ps = 12;
    elseif (i==3)
        SR = 3;
        LR = 15;
        Ps = 6;
    elseif (i==4)
        SR = 3;
        LR = 15;
        Ps = 12;
    end
    for j=1:Ntrials
        if (j < EoR_NL_optimal(i))
            eor{i}(j) = EoR_star(j,j,SR,LR,Ps);
        else
            eor{i}(j) = EoR_star(j,EoR_NL_optimal(i),SR,LR,Ps);
        end
        if (j < f_NL(RoE_NL_optimal(i)))
            roe{i}(j) = RoE_rstar(sum(2:(j+1)),sum(2:(j+1)),SR,LR,Ps)/sum(2:(j+1));
        else
            npresses = RoE_NL_optimal(i) + (j-floor(f_NL(RoE_NL_optimal(i))))*Ps;
            roe{i}(j) = RoE_rstar(npresses,RoE_NL_optimal(i),SR,LR,Ps)/npresses;
        end
    end
end
end

