clear all;
agentOpts.actionSelectionMethod = 'softmax';
agentOpts.updateMethod = 'q-learning';
agentOpts.alpha = .1495;
agentOpts.gamma = .0009;
%agentOpts.epsilon = .01;
agentOpts.maxTemp = .367;
agentOpts.tempDecayFactor = .00;
animalData = load(['processed_data/AZ04_ReProcessedData.mat']); animalData=animalData.ProcessedData;
alphaDecay = .03;
envOpts.itiCost = -1;
envOpts.leverPressCost = 0;
envOpts.rewardFunc = @(u,r) r;
maxPressRemaining = 100;
maxPressAlready = 100;
totalN1 = zeros(maxPressRemaining,maxPressAlready);
totalR1 = zeros(maxPressRemaining,maxPressAlready);
totalN2 = zeros(maxPressRemaining,maxPressAlready);
totalR2 = zeros(maxPressRemaining,maxPressAlready);
totalN3 = zeros(maxPressRemaining,maxPressAlready);
totalR3 = zeros(maxPressRemaining,maxPressAlready);
totalN4 = zeros(maxPressRemaining,maxPressAlready);
totalR4 = zeros(maxPressRemaining,maxPressAlready);
for nAgents=1:10
    agent = leverAgent(agentOpts);
    cumulativeReward = 0;
    for i=1:length(animalData)
        disp(['Day ' num2str(i)])
        envOpts.nPressesForReward_S = animalData{i}.NumPressRequired_S;
        envOpts.nPressesForReward_L = animalData{i}.NumPressRequired_L(1);
        envOpts.Sreward = animalData{i}.Reward_S;
        envOpts.Lreward = animalData{i}.Reward_L;
        env = leverEnv(envOpts,agent);
        ntrials = 0;
        gotReward = [];
        pressInfo = struct('nPressesAlready',{},'nPressesRemaining',{});
        agent.state = 1;
        agent.alpha = 1;
        trialInd = 1;
        step = 1;
        if (strcmp(agent.actionSelectionMethod,'softmax'))
            agent.totalTime = 1;
        end
        while (ntrials <= animalData{i}.TotalTrialsCompleted)
            %disp(['Day ' num2str(i) '    state: ' num2str(agent.state) ' trialInd: ' num2str(trialInd) ' nL: ' num2str(env.nPressesForReward_L)])
            s1 = agent.state;
            states(step) = s1;
            if (s1 == 1)
                ntrials = ntrials+1;
                trialInd = 1;
                nPressesAlready = [];
                nPressesRemaining = [];
            end
            action = agent.act();
            [r,newS] = env.processAction(action);
            if (s1 == 1)
                nPressesAlready(trialInd) = 1;
                if (newS == 2)
                    nPressesRemaining(trialInd) = env.nPressesForReward_S - 1;
                else
                    nPressesRemaining(trialInd) = env.nPressesForReward_L - 1;
                end
                trialInd = trialInd+1;
            else
                if (newS == 1 && action == 2)
                    gotReward(ntrials) = 0;
                    pressInfo(ntrials).nPressesAlready = nPressesAlready;
                    pressInfo(ntrials).nPressesRemaining = nPressesRemaining;
                    %disp(['abort r = ' num2str(r)])
                elseif (newS == 1 && action == 1)
                    gotReward(ntrials) = 1;
                    pressInfo(ntrials).nPressesAlready = nPressesAlready;
                    pressInfo(ntrials).nPressesRemaining = nPressesRemaining;
                    %disp(['success r = ' num2str(r)])
                elseif (newS ~= 1 && action == 1)
                    nPressesAlready(trialInd) = nPressesAlready(trialInd-1)+1;
                    if (newS == 2)
                        nPressesRemaining(trialInd) = env.nPressesForReward_S - nPressesAlready(trialInd);
                    else
                        nPressesRemaining(trialInd) = env.nPressesForReward_L - nPressesAlready(trialInd);
                    end
                    trialInd = trialInd+1;
                    %disp(['lever press  r = ' num2str(r)])
                end
            end
            r = env.rewardFunc(1,r);
            cumulativeReward = cumulativeReward + r;
            agent.state = newS;
            agent.alpha = max(.01,exp(-alphaDecay*ntrials));
            Q(:,:,step) = agent.Q;
            step=step+1;
            agent.updateQ(s1,action,r,newS);
        end
        dayInfo(i).gotReward = gotReward;
        dayInfo(i).pressInfo = pressInfo;
        Qs{i} = Q;
        allStates{i}=states;
        s = sprintf('%1$ixFR%2$i',animalData{i}.Reward_L/animalData{i}.Reward_S,animalData{i}.NumPressRequired_S);
        dayInfo(i).expType = s;
        clear Q;
        clear states;
    end
    disp(['Agent #' num2str(nAgents) '  Cumulative reward = ' num2str(cumulativeReward)])
    for i=1:length(dayInfo)
        N = zeros(maxPressRemaining,maxPressAlready);
        R = zeros(maxPressRemaining,maxPressAlready);
        for j=1:length(dayInfo(i).pressInfo)
            for k=1:length(dayInfo(i).pressInfo(j).nPressesRemaining)
                npr = dayInfo(i).pressInfo(j).nPressesRemaining(k);
                npa = dayInfo(i).pressInfo(j).nPressesAlready(k);
                if (npr > maxPressRemaining || npa > maxPressAlready || npr == 0)
                    continue;
                else
                    N(npr,npa) = N(npr,npa) + 1;
                    R(npr,npa) = R(npr,npa) + dayInfo(i).gotReward(j);
                end
            end
        end
        if (strcmp(dayInfo(i).expType,'2xFR6'))
            totalN1 = totalN1 + N;
            totalR1 = totalR1 + R;
        elseif (strcmp(dayInfo(i).expType,'2xFR12'))
            totalN2 = totalN2 + N;
            totalR2 = totalR2 + R;
        elseif (strcmp(dayInfo(i).expType,'5xFR6'))
            totalN3 = totalN3 + N;
            totalR3 = totalR3 + R;
        elseif (strcmp(dayInfo(i).expType,'5xFR12'))
            totalN4 = totalN4 + N;
            totalR4 = totalR4 + R;
        else
            error(['expType not recognized for day ' num2str(i)])
        end
    end
    clear dayInfo
end
P1 = totalR1./totalN1;
P2 = totalR2./totalN2;
P3 = totalR3./totalN3;
P4 = totalR4./totalN4;
figure;
subplot(2,2,1)
imagesc(P1); set(gca,'ydir','normal'); xlabel('# presses already'); ylabel('# presses remaining'); colormap(jet); caxis([0 1])
title('2xFR6')
subplot(2,2,2)
imagesc(P2); set(gca,'ydir','normal'); xlabel('# presses already'); ylabel('# presses remaining'); colormap(jet); caxis([0 1])
title('2xFR12')
subplot(2,2,3)
imagesc(P3); set(gca,'ydir','normal'); xlabel('# presses already'); ylabel('# presses remaining'); colormap(jet); caxis([0 1])
title('5xFR6')
subplot(2,2,4)
imagesc(P4); set(gca,'ydir','normal'); xlabel('# presses already'); ylabel('# presses remaining'); colormap(jet); caxis([0 1])
title('5xFR12')

titles = {'2xFR6','2xFR12','5xFR6','5xFR12'};
figure;
for i=1:4
    subplot(2,2,i)
    if (i == 1)
        P = P1;
    elseif (i==2)
        P = P2;
    elseif (i==3)
        P = P3;
    elseif (i==4)
        P = P4;
    end
    
    for j=1:10:50
        y = P(:,j); y = y(~isnan(y));
        x = 1:length(y);
        y=y';
        fits{j} = polyfit(x,y,1);
        xs{j} = x;
        ys{j} = y;
    end
    for j=1:length(fits)
        if (length(xs{j}) == 0)
            continue;
        end
        plot(xs{j},xs{j}*fits{j}(1) + fits{j}(2))
        hold on;
    end
    xlabel('# of presses remaining')
    ylabel('P(reward)')
    title(titles{i})
end