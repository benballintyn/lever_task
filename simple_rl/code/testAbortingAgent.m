Q = rand(1000,2);
QSR = rand(12,2);
QPR = rand(1000,2);
alpha1 = .99;
alpha2 = .99;
gamma1 = 0;
gamma2 = .99;
epsilon = .1;

nsessions = 1000;
ntrials = 800;

leverPressCost = 0;
abortCost = 1;

for i=1:nsessions
    sessionType = ceil(rand*4);
    switch sessionType
        case 1
            SR = 3; PR = 6; Ps = 6;
        case 2
            SR = 3; PR = 6; Ps = 12;
        case 3
            SR = 3; PR = 15; Ps = 6;
        case 4
            SR = 3; PR = 15; Ps = 12;
    end
    Pl = 2;
    state1 = 1;
    state2 = 1;
    for j=1:ntrials
        if (rand < epsilon)
            chosenSide = ceil(rand*2);
        else
            [~,chosenSide] = max(Q1(state1,:));
        end
        if (chosenSide == 1) % SR
            completed = 0;
            npress = 0;
            while (~completed)
                if (rand < epsilon)
                    action = ceil(rand*2);
                else
                    [~,action] = max(QSR(state2,:));
                end
                if (action == 1)
                    npress=npress+1;
                    newState2 = state2 + 1;
                    if (npress < Ps)
                        QSR(state2,action) = QSR(state2,action) + alpha2*(leverPressCost + gamma2*max(QSR(newState2,:)) - QSR(state2,action));
                    else
                        QSR(state2,action) = QSR(state2,action) + alpha2*(SR/Ps + gamma2*max(Q(state1,:)) - QSR(state2,action));
                        Q(state1,chosenSide) = Q(state1,chosenSide) + alpha1*(SR/Ps + gamma1*max(Q(state1,:)) - Q(state1,chosenSide));
                        completed = 1;
                    end
                else
                    QSR(state2,action) = QSR(state2,action) + alpha2*(abortCost + gamma2*max(Q(state1,:)) - QSR(state2,action));
                    Q(state1,chosenSide) = Q(state1,chosenSide) + alpha1*(abortCost + gamma1*max(Q(state1,:)) - Q(state1,chosenSide));
                    completed = 1;
                end
                    
            end
        else
            completed = 0;
            npress = 0;
            while (~completed)
                if (rand < epsilon)
                    action = ceil(rand*2);
                else
                    [~,action] = max(QPR(state2,:));
                end
                if (action == 1)
                    npress = npress+1;
                    newState2 = state2 + 1;
                    if (nPress < Pl)
                        QPR(state2,action) = QPR(state2,action) + alpha2*(leverPressCost + gamma2*max(QPR(newState2,:)) - QPR(state2,action));
                    else
                        newState1 = state1 + 1;
                        QPR(state2,action) = QPR(state2,action) + alpha2*(PR/Pl + gamma2*max(Q(newState1,:)) - QPR(state2,action));
                        Q(state1,chosenSide) = Q(state1,chosenSide) + alpha1*(PR/Pl + gamma1*max(Q(newState1,:)) - Q(state1,chosenSide));
                        completed = 1;
                        Pl = Pl + 1;
                    end
                else
                    QPR(state2,action) = QPR(state2,action) + alpha2*(abortCost + gamma2*max(Q(state1,:)) - QPR(state2,action));
                    Q(state1,chosenSide) = Q(state1,chosenSide) + alpha1*(abortCost + gamma1*max(Q(state1,:)) - Q(state1,chosenSide));
                    completed = 1;
                end
            end
        end
            
    end
    disp(['Done with session #' num2str(i)])
end