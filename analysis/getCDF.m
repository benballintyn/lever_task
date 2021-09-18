function [valcdf,valRange] = getCDF(vals,minVal,maxVal,interval)
cdflength = (maxVal - minVal + 1)/interval;
valRange = minVal:interval:maxVal;
valcdf = zeros(1,cdflength);
nvals = length(vals);
for i=1:cdflength
    valcdf(i) = length(find(vals < valRange(i)))/nvals;
end
end

