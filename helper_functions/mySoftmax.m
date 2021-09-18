function [softVals,action] = mySoftmax(vals,temp)
vals = vals-max(vals);
softVals = exp(vals./temp)./sum(exp(vals./temp));
r=rand;
cumVals=cumsum(softVals);
q=find(cumVals > r);
if (isempty(q))
    disp(cumVals)
end
action=q(1);
end

