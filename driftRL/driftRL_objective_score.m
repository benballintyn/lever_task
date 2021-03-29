function [score] = driftRL_objective_score(datadir)
NL_observed = load('optimality/WT/NL_observed.mat'); NL_observed=NL_observed.NL_observed;
nTrialsAborted = load('optimality/WT/nTrialsAborted.mat'); nTrialsAborted=nTrialsAborted.nTrialsAborted;
EoR_optimalities = load('optimality/WT/EoR_optimalities.mat'); EoR_optimalities=EoR_optimalities.EoR_optimalities;
%percentCompletedPR = load('driftRL/abort_analysis/percentCompletedPR.mat'); percentCompletedPR=percentCompletedPR.percentCompletedPR;
percentCompletedPR_allTrials = load('driftRL/abort_analysis/percentCompletedPR_allTrials.mat'); percentCompletedPR_allTrials = percentCompletedPR_allTrials.percentCompletedPR_allTrials;

numLR = load([datadir '/numLR.mat']); numLR=numLR.numLR;
numAborted = load([datadir '/numAborted.mat']); numAborted=numAborted.numAborted;
performanceEoR = load([datadir '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;
agentPercentCompletedPR = load([datadir '/percentCompletedPR.mat']); agentPercentCompletedPR=agentPercentCompletedPR.percentCompletedPR;
%{
for i=1:4
    [mousecdf,valRange] = getCDF(NL_observed{i},1,1000,1); % create cdf of # of NL trials for the mice
    [agentcdf,valRange] = getCDF(numLR{i},1,1000,1);
    wasserstein_distances(i) = wasserstein_1d(mousecdf,agentcdf);
end
%}

logprobs_numLR = zeros(1,4);
logprobs_aborted = zeros(1,4);
logprobs_eor = zeros(1,4);
logprobs_percentCompletedPR = zeros(1,4);
for i=1:4
    [curlogprob_numLR] = logprob2d(NL_observed{i},numLR{i});
    logprobs_numLR(i) = curlogprob_numLR;
    
    [curlogprob_aborted] = logprob2d(nTrialsAborted{i},numAborted{i});
    logprobs_aborted(i) = curlogprob_aborted;
    
    [curlogprob_eor] = logprob2d(EoR_optimalities{i},performanceEoR{i});
    logprobs_eor(i) = curlogprob_eor;
    
    [curlogprob_percentCompletedPR] = logprob2d(percentCompletedPR_allTrials{i},agentPercentCompletedPR{i});
    logprobs_percentCompletedPR(i) = curlogprob_percentCompletedPR;
end

score = -sum([logprobs_numLR logprobs_aborted logprobs_eor logprobs_percentCompletedPR]);
end

