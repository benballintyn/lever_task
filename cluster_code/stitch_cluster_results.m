function [] = stitch_cluster_results(datadir,varargin)
if (~(mod(length(varargin),2) == 0))
    error('Uneven number of varargin fields. For every variable name there must also be a variable range')
end
count = 0;
varCellSizes = [];
for i=1:2:length(varargin)
    count = count+1;
    var_names{count} = varargin{i};
    var_ranges{count} = varargin{i+1};
    varCellSizes = [varCellSizes length(var_ranges{count})];
    params.(var_names{count}) = var_ranges{count};
end
params.var_names = var_names;
fullCellSize = [varCellSizes 4];
folders = dir(datadir);
folders = folders(~ismember({folders.name},{'.','..'}));
folders = folders([folders.isdir]);

performanceEoR = cell(fullCellSize); % Fraction of max EoR values
performanceRoE = cell(fullCellSize);
nL = cell(fullCellSize); % # of choices of PR side for each timestep
nS = cell(fullCellSize); % # of choices of SR side for each timestep
pL = cell(fullCellSize); % P(PR) for each timestep. computed from nL and nS
pS = cell(fullCellSize); % P(SR) for each timestep. computed from nL and nS
numLR = cell(fullCellSize); % Total # of PR trials
numSR = cell(fullCellSize); % Total # of SR trials
startTime = tic;
nFolders = length(folders);
for i=1:nFolders
    curRunParams = load([datadir '/' folders(i).name '/runParams.mat']); curRunParams=curRunParams.runParams;
    if (i == 1)
        sweepParams.agentType = curRunParams.agentType;
        sweepParams.actionSelectionMethod = curRunParams.actionSelectionMethod;
        sweepParams.utilityFunc1 = curRunParams.utilityFunc1;
        sweepParams.utilityFunc2 = curRunParams.utilityFunc2;
    else
        if (~strcmp(sweepParams.agentType,curRunParams.agentType))
            error(['agentType from ' folders(i).name ' does not match that from ' folders(1).name])
        end
        if (~strcmp(sweepParams.actionSelectionMethod,curRunParams.actionSelectionMethod))
            error(['actionSelectionMethod from ' folders(i).name ' does not match that from ' folders(1).name])
        end
        if (~strcmp(sweepParams.utilityFunc1,curRunParams.utilityFunc1))
            error(['utilityFunc1 from ' folders(i).name ' does not match that from ' folders(1).name])
        end
        if (~strcmp(sweepParams.utilityFunc2,curRunParams.utilityFunc2))
            error(['utilityFunc2 from ' folders(i).name ' does not match that from ' folders(1).name])
        end
    end
    agentParams = load([datadir '/' folders(i).name '/agentParams.mat']); agentParams=agentParams.agentParams;
    var_val = zeros(1,length(var_names));
    cellInds = zeros(1,length(var_names));
    for j=1:length(var_names)
        var_val(j) = agentParams.(var_names{j});
        cellInds(j) = find(var_ranges{j} == var_val(j));
    end
    tmpperformanceEoR = load([datadir '/' folders(i).name '/performanceEoR.mat']); tmpperformanceEoR=tmpperformanceEoR.performanceEoR;
    tmpperformanceRoE = load([datadir '/' folders(i).name '/performanceRoE.mat']); tmpperformanceRoE=tmpperformanceRoE.performanceRoE;
    tmpnL = load([datadir '/' folders(i).name '/nL.mat']); tmpnL=tmpnL.nL;
    tmpnS = load([datadir '/' folders(i).name '/nS.mat']); tmpnS=tmpnS.nS;
    tmppL = load([datadir '/' folders(i).name '/pL.mat']); tmppL=tmppL.pL;
    tmppS = load([datadir '/' folders(i).name '/pS.mat']); tmppS=tmppS.pS;
    tmpnumLR = load([datadir '/' folders(i).name '/numLR.mat']); tmpnumLR=tmpnumLR.numLR;
    tmpnumSR = load([datadir '/' folders(i).name '/numSR.mat']); tmpnumSR=tmpnumSR.numSR;
    for j=1:4
        linearIndex = sub2ind_cell(nL,[cellInds j]);
        performanceEoR{linearIndex} = tmpperformanceEoR{j};
        performanceRoE{linearIndex} = tmpperformanceRoE{j};
        nL{linearIndex} = tmpnL{j};
        nS{linearIndex} = tmpnS{j};
        pL{linearIndex} = tmppL{j};
        pS{linearIndex} = tmppS{j};
        numLR{linearIndex} = tmpnumLR{j};
        numSR{linearIndex} = tmpnumSR{j};
    end
    if (mod(i,1000) == 0)
        timeElapsed = toc(startTime);
        speed = timeElapsed/i;
        timeRemaining_mins = speed*(nFolders - i)/60;
        disp(['Done with first ' num2str(i) ' folders. ~' num2str(timeRemaining_mins) ' minutes remaining'])
    end
end
save([datadir '/performanceEoR.mat'],'performanceEoR','-mat')
save([datadir '/performanceRoE.mat'],'performanceRoE','-mat')
save([datadir '/nL.mat'],'nL','-mat')
save([datadir '/nS.mat'],'nS','-mat')
save([datadir '/pL.mat'],'pL','-mat')
save([datadir '/pS.mat'],'pS','-mat')
save([datadir '/numLR.mat'],'numLR','-mat')
save([datadir '/numSR.mat'],'numSR','-mat')
save([datadir '/params.mat'],'params','-mat')
save([datadir '/sweepParams.mat'],'sweepParams','-mat')
end

