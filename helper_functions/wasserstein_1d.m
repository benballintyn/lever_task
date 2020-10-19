function [wd] = wasserstein_1d(cdf1,cdf2)
wd = sum(abs(cdf1-cdf2));
end

