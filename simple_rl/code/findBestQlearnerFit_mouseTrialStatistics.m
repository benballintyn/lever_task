clear all;
% compare rl nL distributions to mice nL distributions using Wasserstein
% Distance
mice = {'AZ04','AZ08','HP01','HP02','HP03','HP04','MA01','NS07','NS09','NS10'};
NL_observed = load('optimality/WT/NL_observed.mat'); NL_observed=NL_observed.NL_observed;
datadir = '~/phd/lever_task/simple_rl/utility_func1_softmax_largeSweep_mouseTrialStatistics';
nL = load([datadir '/nL.mat']); nL=nL.nL;
numLR = load([datadir '/numLR.mat']); numLR=numLR.numLR;
sessionTypes = {'2xFR6','2xFR12','5xFR6','5xFR12'};
use_e_greedy = 0;
alphas = load([datadir '/alphas.mat']); alphas=alphas.alphas;
gammas = load([datadir '/gammas.mat']); gammas=gammas.gammas;
performanceEoR = load([datadir '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;
performanceRoE = load([datadir '/performanceRoE.mat']); performanceRoE=performanceRoE.performanceRoE;
if (use_e_greedy)
    epsilons = load([datadir '/epsilons.mat']); epsilons=epsilons.epsilons;
else
    temps = load([datadir '/temps.mat']); temps=temps.temps;
end

n1 = size(numLR,1); % # of alpha values
n2 = size(numLR,2); % # of gamma values
n3 = size(numLR,3); % # of epsilon/temmp values
n4 = size(numLR,4); % # of session types
for i=1:4
    [valcdf,valRange] = getCDF(NL_observed{i},1,1000,1); % create cdf of # of NL trials for the mice
    mouseCDF{i} = valcdf;
end
% For each paramter set and session type, compute the RL agent's cdf and
% compute the Wasserstein distance with the mouse distribution from the
% corresponding session type
rlcdfs = cell(n1,n2,n3,n4); 
for i=1:n1
    for j=1:n2
        for k=1:n3
            for l=1:n4
                [valcdf,valRange] = getCDF(numLR{i,j,k,l},1,1000,1);
                rlcdfs{i,j,k,l} = valcdf;
                wd(i,j,k,l) = wasserstein_1d(mouseCDF{l},valcdf);
            end
        end
    end
end
% Go through all of the wasserstein distances and find the parameter set
% with the lowest. Also use the RoE optimality values of each rl agent and
% find the parameter set that his the highest RoE optimality.
bestInds = zeros(4,3);
bestPerformanceInds = zeros(4,3);
for i=1:n4
    minwd = double(intmax);
    maxEoR = 0;
    for j=1:n1
        for k=1:n2
            for l=1:n3
                if (wd(j,k,l,i) < minwd)
                    minwd = wd(j,k,l,i);
                    bestInds(i,:) = [j k l];
                end
                if (mean(performanceEoR{j,k,l,i}) > maxEoR)
                    maxEoR = mean(performanceEoR{j,k,l,i});
                    bestPerformanceInds(i,:) = [j k l];
                end
            end
        end
    end
end
save([datadir '/bestInds.mat'],'bestInds','-mat')
save([datadir '/bestPerformanceInds.mat'],'bestPerformanceInds','-mat')

if (~exist([datadir '/figures'],'dir'))
    mkdir([datadir '/figures'])
end

% Plot 1 subplot per session type showing the #NL trial distribution for
% mice as well as the RL agent that has the best fit
figure;
for i=1:4
    subplot(2,2,i)
    histogram(NL_observed{i},'binwidth',2,'normalization','pdf');
    hold on;
    histogram(numLR{bestInds(i,1),bestInds(i,2),bestInds(i,3),i},'binwidth',2,'normalization','pdf')
    if (use_e_greedy)
        title([sessionTypes{i} ' | \alpha = ' num2str(alphas(bestInds(i,1))) ' \gamma = ' num2str(gammas(bestInds(i,2))) '\epsilon = ' num2str(epsilons(bestInds(i,3)))])
    else
        title([sessionTypes{i} ' | \alpha = ' num2str(alphas(bestInds(i,1))) ' \gamma = ' num2str(gammas(bestInds(i,2))) '\tau = ' num2str(temps(bestInds(i,3)))])
    end
end
suptitle('Best match to mice')
set(gcf,'Position',[10 10 1400 1400])
saveas(gcf,[datadir '/figures/num_PR_trials_best_fit.fig'],'fig')
saveas(gcf,[datadir '/figures/num_PR_trials_best_fit.eps'],'eps')


figure;
for i=1:4
    subplot(2,2,i)
    histogram(NL_observed{i},'binwidth',2,'normalization','pdf');
    hold on;
    histogram(numLR{bestPerformanceInds(i,1),bestPerformanceInds(i,2),bestPerformanceInds(i,3),i},'binwidth',2,'normalization','pdf')
    if (use_e_greedy)
        title([sessionTypes{i} ' | \alpha = ' num2str(alphas(bestPerformanceInds(i,1))) ' \gamma = ' num2str(gammas(bestPerformanceInds(i,2))) '\epsilon = ' num2str(epsilons(bestPerformanceInds(i,3)))])
    else
        title([sessionTypes{i} ' | \alpha = ' num2str(alphas(bestPerformanceInds(i,1))) ' \gamma = ' num2str(gammas(bestPerformanceInds(i,2))) '\tau = ' num2str(temps(bestPerformanceInds(i,3)))])
    end
end
suptitle('Best overall agent')
set(gcf,'Position',[10 10 1400 1400])
saveas(gcf,[datadir '/figures/num_PR_trials_best_performance.fig'],'fig')
saveas(gcf,[datadir '/figures/num_PR_trials_best_performance.eps'],'eps')

figure;
for i=1:4
    subplot(2,2,i)
    surf(gammas(1:end-1),alphas,squeeze(wd(:,1:end-1,bestInds(i,3),i)))
    xlabel('\gamma')
    ylabel('\alpha')
    if (use_e_greedy)
        title(['\epsilon = ' num2str(epsilons(bestInds(i,3)))])
    else
        title(['\tau = ' num2str(temps(bestInds(i,3)))])
    end
end
set(gcf,'Position',[10 10 1400 1400])
saveas(gcf,[datadir '/figures/alpha_gamma_slice.fig'],'fig')
saveas(gcf,[datadir '/figures/alpha_gamma_slice.eps'],'eps')

figure;
for i=1:4
    subplot(2,2,i)
    if (use_e_greedy)
        surf(epsilons,alphas,squeeze(wd(:,bestInds(i,2),:,i)))
        xlabel('\epsilon')
    else
        surf(temps,alphas,squeeze(wd(:,bestInds(i,2),:,i)))
        xlabel('\tau')
    end
    ylabel('\alpha')
    title(['\gamma = ' num2str(gammas(bestInds(i,2)))])
end
set(gcf,'Position',[10 10 1400 1400])
if (use_e_greedy)
    saveas(gcf,[datadir '/figures/alpha_epsilon_slice.fig'],'fig')
    saveas(gcf,[datadir '/figures/alpha_epsilon_slice.eps'],'eps')
else
    saveas(gcf,[datadir '/figures/alpha_temp_slice.fig'],'fig')
    saveas(gcf,[datadir '/figures/alpha_temp_slice.eps'],'eps')
end

figure;
for i=1:4
    subplot(2,2,i)
    if (use_e_greedy)
        surf(epsilons,gammas(1:end-1),squeeze(wd(bestInds(i,1),1:end-1,:,i)))
        xlabel('\epsilon')
    else
        surf(temps,gammas(1:end-1),squeeze(wd(bestInds(i,1),1:end-1,:,i)))
        xlabel('\tau')
    end
    title(['\alpha = ' num2str(alphas(bestInds(i,1)))])
    ylabel('\gamma')
end
set(gcf,'Position',[10 10 1400 1400])
if (use_e_greedy)
    saveas(gcf,[datadir '/figures/gamma_epsilon_slice.fig'],'fig')
    saveas(gcf,[datadir '/figures/gamma_epsilon_slice.eps'],'eps')
else
    saveas(gcf,[datadir '/figures/gamma_temp_slice.fig'],'fig')
    saveas(gcf,[datadir '/figures/gamma_temp_slice.eps'],'eps')
end
