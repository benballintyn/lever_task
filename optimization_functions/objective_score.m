function [score] = objective_score(allActionProbs)
f = @(actionProbs,A) 1./(1 + exp(-A*tanh((1/A)*log(actionProbs./(1 - actionProbs)))));
score = sum(-log(f(allActionProbs,10)));
end

