function [logprob,pd,xi] = logprob2d(actualDistribution,modelDistribution,varargin)
p=inputParser;
addRequired(p,'actualDistribution')
addRequired(p,'modelDistribution')
addParameter(p,'bandwidth',nan,@isnumeric)
parse(p,actualDistribution,modelDistribution,varargin{:});
if (isnan(p.Results.bandwidth))
    [pd,xi] = ksdensity(modelDistribution,actualDistribution,'NumPoints',10000);
else
    [pd,xi] = ksdensity(modelDistribution,actualDistribution,'NumPoints',10000,'bandwidth',p.Results.bandwidth);
end
pd(pd == 0) = 1/realmax;
logprob = sum(log(pd));
end

