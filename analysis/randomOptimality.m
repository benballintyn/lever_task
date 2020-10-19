function [RoEoptimalities,EoRoptimalities] = randomOptimality(maxTrials,nrepetitions)
RoE_rstar = @(N,NL,SR,LR,Ps) (N - NL)*(SR/Ps) + LR*(.5*(sqrt(8*NL+9)-3));
RoE_NL_optimal = [70.875 286.875 448.875 1798.875];
EoR_NL_optimal = [10 22 28 58];
eg = double(vpa(eulergamma));
EoR_star = @(N,NL,SR,LR,Ps) (1/N)*(LR*(psi(0,NL+2) + eg - 1) + (N - NL)*(SR/Ps));

RoEoptimalities = cell(1,4);
EoRoptimalities = cell(1,4);
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
    RoEoptimalities{i} = zeros(maxTrials,nrepetitions);
    EoRoptimalities{i} = zeros(maxTrials,nrepetitions);
    for j=1:maxTrials
        disp(['Done with sessionType ' num2str(i) ' numTrials = ' num2str(j)])
        parfor k=1:nrepetitions
            rewards = zeros(1,j);
            presses = zeros(1,j);
            Pl = 2;
            for l=1:j
                r = rand;
                if (r < .5)
                    rewards(l) = SR;
                    presses(l) = Ps;
                else
                    rewards(l) = LR;
                    presses(l) = Pl;
                    Pl=Pl+1;
                end
            end
            if (sum(presses) < RoE_NL_optimal(i))
                bestNL = sum(presses);
            else
                bestNL = RoE_NL_optimal(i);
            end
            if (j < EoR_NL_optimal(i))
                curEoRstar = EoR_star(j,j,SR,LR,Ps);
            else
                curEoRstar = EoR_star(j,EoR_NL_optimal(i),SR,LR,Ps);
            end
            curRstar = RoE_rstar(sum(presses),bestNL,SR,LR,Ps);
            roes(k) = sum(rewards)/curRstar;
            EoRobserved = mean(rewards./presses);
            eors(k) = EoRobserved/curEoRstar;
            %randOptimalities{i}(j,k) = sum(rewards)/curRstar;
        end
        RoEoptimalities{i}(j,:) = roes;
        EoRoptimalities{i}(j,:) = eors;
    end
end
end

