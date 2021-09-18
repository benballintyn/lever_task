function [p,H,p2] = scheirer_ray_hare(data,group)

[R] = tiedrank(data);

% compute the anova SS's on the ranks:
[~,p2] = anovan(R,group,'display','off','model','interaction');
nvar = size(p2,1)-2;
for i=1:nvar
    df(i) = p2{i+1,3};
    SS(i) = p2{i+1,2};
end
% total MS
MS = sum(SS)/sum(df);
% compute ratios
for i=1:nvar-1
    H(i) = SS(i)/MS;
end

% these are compared to a chi-2-distribution with df deg of freedom:
for i=1:nvar-1
    p(i) = 1-chi2cdf(H(i),df(i));
end