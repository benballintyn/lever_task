clear;

% datapath = 'C:\Users\nategr8\Documents\Incentive Motivation\data\NewParam\WT_data\';
datapath = '~/phd/lever_task/raw_data/';

ntrials_PR_0606 = zeros(10,100);
ntrials_PR_1206 = zeros(10,100);
ntrials_PR_0615 = zeros(10,100);
ntrials_PR_1215 = zeros(10,100);

comp_PR_0606 = zeros(10,100);
comp_PR_1206 = zeros(10,100);
comp_PR_0615 = zeros(10,100);
comp_PR_1215 = zeros(10,100);


flist = dir([datapath '*.mat']);
for f = 1:length(flist)
    load([datapath flist(f).name]);

    ntrials = SessionData.nTrials;
    ntrials_PR = zeros(10,100);
    FR_press_req = SessionData.TrialSettings(1).GUI.NumPressRequired_S;
    PR_press_req_max = SessionData.TrialSettings(end).GUI.NumPressRequired_L;
    PR_volume = SessionData.TrialSettings(1).GUI.Reward_L;
    comp_PR = zeros(10,100);


    for i = 1:ntrials
        PR_press_req = SessionData.TrialSettings(i).GUI.NumPressRequired_L;
        states = SessionData.RawEvents.Trial{1,i}.States;
        if ~isnan(SessionData.RawEvents.Trial{1,i}.States.Side_L(1))
            ntrials_PR(1,PR_press_req-1) = ntrials_PR(1,PR_press_req-1) + 1;
            if ~isnan(states.WaitForPoke_L())
                comp_PR(1,PR_press_req-1) = comp_PR(1,PR_press_req-1) + 1;
            end
            if PR_press_req > 5
                if ~isnan(states.WaitForPress_L6(1))
                    ntrials_PR(2,PR_press_req - 5) = ntrials_PR(2,PR_press_req - 5) + 1;
                    if ~isnan(states.WaitForPoke_L())
                        comp_PR(2,PR_press_req - 5) = comp_PR(2,PR_press_req - 5) + 1;
                    end
                end
            end
            if PR_press_req > 10
                if ~isnan(states.WaitForPress_L11(1))
                    ntrials_PR(3,PR_press_req - 10) = ntrials_PR(3,PR_press_req - 10) + 1;
                    if ~isnan(states.WaitForPoke_L())
                        comp_PR(3,PR_press_req - 10) = comp_PR(3,PR_press_req - 10) + 1;
                    end
                end  
            end
            if PR_press_req > 15
                if ~isnan(states.WaitForPress_L16(1))
                    ntrials_PR(4,PR_press_req - 15) = ntrials_PR(4,PR_press_req - 15) + 1;
                    if ~isnan(states.WaitForPoke_L())
                        comp_PR(4,PR_press_req - 15) = comp_PR(4,PR_press_req - 15) + 1;
                    end
                end
            end
            if PR_press_req > 20
                if ~isnan(states.WaitForPress_L21(1))
                    ntrials_PR(5,PR_press_req - 20) = ntrials_PR(5,PR_press_req - 20) + 1;
                    if ~isnan(states.WaitForPoke_L())
                        comp_PR(5,PR_press_req - 20) = comp_PR(5,PR_press_req - 20) + 1;
                    end
                end
            end
            if PR_press_req > 25
                if ~isnan(states.WaitForPress_L26(1))
                    ntrials_PR(6,PR_press_req - 25) = ntrials_PR(6,PR_press_req - 25) + 1;
                    if ~isnan(states.WaitForPoke_L())
                        comp_PR(6,PR_press_req - 25) = comp_PR(6,PR_press_req - 25) + 1;
                    end
                end
            end
            if PR_press_req > 30
                if ~isnan(states.WaitForPress_L31(1))
                    ntrials_PR(7,PR_press_req - 30) = ntrials_PR(7,PR_press_req - 30) + 1;
                    if ~isnan(states.WaitForPoke_L())
                        comp_PR(7,PR_press_req - 30) = comp_PR(7,PR_press_req - 30) + 1;
                    end
                end
            end
            if PR_press_req > 35
                if ~isnan(states.WaitForPress_L36(1))
                    ntrials_PR(8,PR_press_req - 35) = ntrials_PR(8,PR_press_req - 35) + 1;
                    if ~isnan(states.WaitForPoke_L())
                        comp_PR(8,PR_press_req - 35) = comp_PR(8,PR_press_req - 35) + 1;
                    end
                end
            end
            if PR_press_req > 40
                if ~isnan(states.WaitForPress_L41(1))
                    ntrials_PR(9,PR_press_req - 40) = ntrials_PR(9,PR_press_req - 40) + 1;
                    if ~isnan(states.WaitForPoke_L())
                        comp_PR(9,PR_press_req - 40) = comp_PR(9,PR_press_req - 40) + 1;
                    end
                end
            end
            if PR_press_req > 45
                if ~isnan(states.WaitForPress_L46(1))
                    ntrials_PR(10,PR_press_req - 45) = ntrials_PR(10,PR_press_req - 45) + 1;
                    if ~isnan(states.WaitForPoke_L())
                        comp_PR(10,PR_press_req - 45) = comp_PR(10,PR_press_req - 45) + 1;
                    end
                end
            end  
        end
    end

    if FR_press_req == 6 && PR_volume == 6
        ntrials_PR_0606 = ntrials_PR_0606 + ntrials_PR;
        comp_PR_0606 = comp_PR_0606 + comp_PR;
    elseif FR_press_req == 12 && PR_volume == 6
        ntrials_PR_1206 = ntrials_PR_1206 + ntrials_PR;
        comp_PR_1206 = comp_PR_1206 + comp_PR;
    elseif FR_press_req == 6 && PR_volume == 15
        ntrials_PR_0615 = ntrials_PR_0615 + ntrials_PR;
        comp_PR_0615 = comp_PR_0615 + comp_PR;
    elseif FR_press_req == 12 && PR_volume == 15
        ntrials_PR_1215 = ntrials_PR_1215 + ntrials_PR;
        comp_PR_1215 = comp_PR_1215 + comp_PR;
    end
end

p_earnPR_0606 = comp_PR_0606./ntrials_PR_0606;
p_earnPR_1206 = comp_PR_1206./ntrials_PR_1206;
p_earnPR_0615 = comp_PR_0615./ntrials_PR_0615;
p_earnPR_1215 = comp_PR_1215./ntrials_PR_1215;

colors = rgb('black','crimson','dark orange','goldenrod','lime green',...
'kelley green','bright blue','dark blue','purple','lilac');

%%
fit0606 = nan(3,100,10,10);
figure; hold on
title('2xFR6')
ylabel('p(earn)')
xlabel('Presses remaining')
for i = 1:10
    ntrials_toplot = 0;
    for j = 1:100
        if ntrials_PR_0606(i,j) >= 10
            ntrials_toplot = ntrials_toplot + 1;
        end
    end
    if ntrials_toplot >= 10
        plot(p_earnPR_0606(i,1:ntrials_toplot), 'o', 'Color', colors(i,:))
        [b0606(:,i), bint0606(:,:,i)] = regress(p_earnPR_0606(i,1:ntrials_toplot)',[[1:ntrials_toplot]' ones(ntrials_toplot,1)]);
        fit0606(1,1:ntrials_toplot,i) = polyval(b0606(:,i),1:ntrials_toplot);
        fit0606(2,1:ntrials_toplot,i) = polyval(bint0606(:,1,i),1:ntrials_toplot);
        fit0606(3,1:ntrials_toplot,i) = polyval(bint0606(:,2,i),1:ntrials_toplot);
        plot(fit0606(1,:,i),'-', 'Color', colors(i,:))
        fill([1:ntrials_toplot flip(1:ntrials_toplot)]',[fit0606(2,1:ntrials_toplot,i) flip(fit0606(3,1:ntrials_toplot,i))]', colors(i,:), 'facealpha', 0.4, 'edgecolor', 'none');
    end
end

figure;
for i = 1:size(b0606,2)
    bar(i,b0606(1,i),'FaceColor',colors(i,:))
    hold on
end
title('2xFR6')
ylabel('Slope of P(earn)')

%% 
fit1206 = nan(3,100,10);
figure; hold on
title('2xFR12')
ylabel('p(earn)')
xlabel('Presses remaining')
for i = 1:10
    ntrials_toplot = 0;
    for j = 1:100
        if ntrials_PR_1206(i,j) >= 10
            ntrials_toplot = ntrials_toplot + 1;
        end
    end
    if ntrials_toplot >= 10
        plot(p_earnPR_1206(i,1:ntrials_toplot), 'o', 'Color', colors(i,:))
        [b1206(:,i), bint1206(:,:,i)] = regress(p_earnPR_1206(i,1:ntrials_toplot)',[[1:ntrials_toplot]' ones(ntrials_toplot,1)]);
        fit1206(1,1:ntrials_toplot,i) = polyval(b1206(:,i),1:ntrials_toplot);
        fit1206(2,1:ntrials_toplot,i) = polyval(bint1206(:,1,i),1:ntrials_toplot);
        fit1206(3,1:ntrials_toplot,i) = polyval(bint1206(:,2,i),1:ntrials_toplot);
        plot(fit1206(1,:,i),'-', 'Color', colors(i,:))
        fill([1:ntrials_toplot flip(1:ntrials_toplot)]',[fit1206(2,1:ntrials_toplot,i) flip(fit1206(3,1:ntrials_toplot,i))]', colors(i,:), 'facealpha', 0.4, 'edgecolor', 'none');
    end
end

figure;
for i = 1:size(b1206,2)
    bar(i,b1206(1,i),'FaceColor',colors(i,:))
    hold on
end
title('2xFR12')
ylabel('Slope of P(earn)')

%% 
fit0615 = nan(3,100,10);
figure; hold on
title('5xFR6')
ylabel('p(earn)')
xlabel('Presses remaining')
for i = 1:10
    ntrials_toplot = 0;
    for j = 1:100
        if ntrials_PR_0615(i,j) >= 10
            ntrials_toplot = ntrials_toplot + 1;
        end
    end
    if ntrials_toplot >= 10
        plot(p_earnPR_0615(i,1:ntrials_toplot), 'o', 'Color', colors(i,:))
        [b0615(:,i), bint0615(:,:,i)] = regress(p_earnPR_0615(i,1:ntrials_toplot)',[[1:ntrials_toplot]' ones(ntrials_toplot,1)]);
        fit0615(1,1:ntrials_toplot,i) = polyval(b0615(:,i),1:ntrials_toplot);
        fit0615(2,1:ntrials_toplot,i) = polyval(bint0615(:,1,i),1:ntrials_toplot);
        fit0615(3,1:ntrials_toplot,i) = polyval(bint0615(:,2,i),1:ntrials_toplot);
        plot(fit0615(1,:,i),'-', 'Color', colors(i,:))
        fill([1:ntrials_toplot flip(1:ntrials_toplot)]',[fit0615(2,1:ntrials_toplot,i) flip(fit0615(3,1:ntrials_toplot,i))]', colors(i,:), 'facealpha', 0.4, 'edgecolor', 'none');
    end
end

figure;
for i = 1:size(b0615,2)
    bar(i,b0615(1,i),'FaceColor',colors(i,:))
    hold on
end
title('5xFR6')
ylabel('Slope of P(earn)')

%% 
fit1215 = nan(3,100,10);
figure; hold on
title('5xFR12')
ylabel('p(earn)')
xlabel('Presses remaining')
for i = 1:10
    ntrials_toplot = 0;
    for j = 1:100
        if ntrials_PR_1215(i,j) >= 10
            ntrials_toplot = ntrials_toplot + 1;
        end
    end
    if ntrials_toplot >= 10
        plot(p_earnPR_1215(i,1:ntrials_toplot), 'o', 'Color', colors(i,:))
        [b1215(:,i), bint1215(:,:,i)] = regress(p_earnPR_1215(i,1:ntrials_toplot)',[[1:ntrials_toplot]' ones(ntrials_toplot,1)]);
        fit1215(1,1:ntrials_toplot,i) = polyval(b1215(:,i),1:ntrials_toplot);
        fit1215(2,1:ntrials_toplot,i) = polyval(bint1215(:,1,i),1:ntrials_toplot);
        fit1215(3,1:ntrials_toplot,i) = polyval(bint1215(:,2,i),1:ntrials_toplot);
        plot(fit1215(1,:,i),'-', 'Color', colors(i,:))
        fill([1:ntrials_toplot flip(1:ntrials_toplot)]',[fit1215(2,1:ntrials_toplot,i) flip(fit1215(3,1:ntrials_toplot,i))]', colors(i,:), 'facealpha', 0.4, 'edgecolor', 'none');
    end
end

figure;
for i = 1:size(b1215,2)
    bar(i,b1215(1,i),'FaceColor',colors(i,:))
    hold on
end
title('5xFR12')
ylabel('Slope of P(earn)')
