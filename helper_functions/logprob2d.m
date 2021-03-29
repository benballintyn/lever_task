function [logprob,pd] = logprob2d(actualDistribution,modelDistribution)
[pd,xi] = ksdensity(modelDistribution,actualDistribution,'NumPoints',10000);
pd(pd == 0) = 1/realmax;
logprob = sum(log(pd));
end

