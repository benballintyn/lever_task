function [counts3D] = make3Dcounts(datadir,varargin)
p=inputParser;
addRequired(p,'datadir',@ischar)
addParameter(p,'nPR_counts',15,@isnumeric)
addParameter(p,'nPR_binSize',5,@isnumeric)
addParameter(p,'nAborted_counts',14,@isnumeric)
addParameter(p,'nAborted_binSize',5,@isnumeric)
addParameter(p,'EoR_counts',20,@isnumeric)
addParameter(p,'EoR_binSize',.025,@isnumeric)
parse(p,datadir,varargin{:})

nPR_counts = p.Results.nPR_counts;
nPR_binSize = p.Results.nPR_binSize;
nAborted_counts = p.Results.nAborted_counts;
nAborted_binSize = p.Results.nAborted_binSize;
EoR_counts = p.Results.EoR_counts;
EoR_binSize = p.Results.EoR_binSize;

numLR = load([datadir '/numLR.mat']); numLR=numLR.numLR;
numAborted = load([datadir '/numAborted.mat']); numAborted=numAborted.numAborted;
performanceEoR = load([datadir '/performanceEoR.mat']); performanceEoR=performanceEoR.performanceEoR;

counts3D = cell(1,4);
for i=1:4
    counts3D{i} = zeros(nPR_counts,nAborted_counts,EoR_counts);
end

for i=1:4
    nDataPoints = length(numLR{i});
    if (length(numAborted{i}) ~= nDataPoints || length(performanceEoR{i}) ~= nDataPoints)
        error(['Inconsistent number of data points for session type ' num2str(i)])
    end
    for j=1:nDataPoints
        % Get index for # of PR trials completed
        PR_count_ind = ceil(numLR{i}(j)/nPR_binSize);
        if (PR_count_ind == 0)
            PR_count_ind = 1;
        elseif (PR_count_ind > nPR_counts)
            PR_count_ind = nPR_counts;
        end
        
        % Get index for # of trials aborted
        nAborted_ind = ceil(numAborted{i}(j)/nAborted_binSize);
        if (nAborted_ind == 0)
            nAborted_ind = 1;
        elseif (nAborted_ind > nAborted_counts)
            nAborted_ind = nAborted_counts;
        end
        
        % Get index for EoR optimality
        EoR_ind = ceil(performanceEoR{i}(j)/EoR_binSize) - 20; % -20 is because we don't use bins for 0-.5
        if (EoR_ind <= 0)
            EoR_ind = 1;
        elseif (EoR_ind > EoR_counts)
            EoR_ind = EoR_counts;
        end
        
        % add data point to joint count matrix
        counts3D{i}(PR_count_ind,nAborted_ind,EoR_ind) = ...
            counts3D{i}(PR_count_ind,nAborted_ind,EoR_ind) + 1;
    end
end
end

