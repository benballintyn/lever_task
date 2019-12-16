function [softVals,action] = mySoftmax(vals,temp)
vals = vals-max(vals);
softVals = exp(vals./temp)./sum(exp(vals./temp));
r=rand;
cumVals=cumsum(softVals);
q=find(cumVals > r);
action=q(1);
end

