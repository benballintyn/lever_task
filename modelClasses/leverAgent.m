classdef leverAgent < matlab.mixin.Copyable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties
        actionSelectionMethod
        updateMethod
        valueAdjustment
        utilityFunc
        state
        alpha
        gamma
        epsilon
        c
        maxTemp
        tempDecayFactor
        tempFunc
        Q
        locationModel
        actionCount
        totalTime
    end
    
    methods
        function obj = leverAgent(options)
            % enforce the presence of certain inputs
            if (~isfield(options,'actionSelectionMethod'))
                error('No actionSelectionMethod specified')
            end
            if (~isfield(options,'updateMethod'))
                error('No update method specified')
            end
            % define the default parameters
            defOpts = struct('actionSelectionMethod','e-greedy', ...
                             'alpha',.01, ...
                             'gamma',.9, ...
                             'epsilon',.05, ...
                             'tempFunc',@(t) max(.001,obj.maxTemp*exp(-obj.tempDecayFactor*t)), ...
                             'utilityFunc',@(t) 1, ...
                             'maxTemp',1, ...
                             'tempDecayFactor',.001, ...
                             'valueAdjustment','standard', ...
                             'c', 1);
            % begin building agent from inputs
            if (strcmp(options.actionSelectionMethod,'e-greedy'))
                obj.actionSelectionMethod = 'e-greedy';
                if (isfield(options,'epsilon'))
                    obj.epsilon = options.epsilon;
                else
                    disp('No epsilon provided... using default value')
                    obj.epsilon = defOpts.epsilon;
                end
            end
            if (strcmp(options.actionSelectionMethod,'softmax'))
                obj.actionSelectionMethod = 'softmax';
                if (isfield(options,'tempFunc'))
                    obj.tempFunc = options.tempFunc;
                else
                    disp('No tempFunc provided... using default value')
                    obj.tempFunc = defOpts.tempFunc;
                end
                if (isfield(options,'maxTemp'))
                    obj.maxTemp = options.maxTemp;
                else
                    disp('No maxTemp provided... using default value')
                    obj.maxTemp = defOpts.maxTemp;
                end
                
                if (isfield(options,'tempDecayFactor'))
                    obj.tempDecayFactor = options.tempDecayFactor;
                else
                    disp('No tempDecayFactor provided... using default value')
                    obj.tempDecayFactor = defOpts.tempDecayFactor;
                end
            end
            
            % set value adjustment method
            if (isfield(options,'valueAdjustment'))
                obj.valueAdjustment = options.valueAdjustment;
            else
                disp('No value adjustment method provided... using default')
                obj.valueAdjustment = defOpts.valueAdjustment;
            end
            if (strcmp(obj.valueAdjustment,'UCB'))
                if (isfield(options,'c'))
                    obj.c = options.c;
                else
                    disp('No c provided... using default value')
                    obj.c = defOpts.valueAdjustment;
                end
            end
            % initialize Q table and action count table
            if (strcmp(options.updateMethod,'q-learning'))
                obj.updateMethod = 'q-learning';
                obj.Q = rand(3,2);
                obj.actionCount = zeros(3,2);
            end
            
            % set utility function (required input)
            if (~isfield(options,'utilityFunc'))
                disp('No utilityFunc provided... using default')
                obj.utilityFunc = defOpts.utilityFunc;
            else
                obj.utilityFunc = options.utilityFunc;
            end
            
            % Initialize location, utility models
            if (isfield(options,'useModel') && options.useModel == 1)
               % not implemented
            end
            
            % initialize alpha and gamma parameters
            if (isfield(options,'alpha'))
                obj.alpha = options.alpha;
            else
                disp('No alpha provided... using default value')
                obj.alpha = defOpts.alpha;
            end
            if (isfield(options,'gamma'))
                obj.gamma = options.gamma;
            else
                disp('No gamma provided... using default value')
                obj.gamma = defOpts.gamma;
            end
            % initialize state
            obj.state = 1;
            
            %initialize totalTime
            obj.totalTime = 1;
        end
        
        function action = act(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.totalTime = obj.totalTime+1;
            s = obj.state;
            if (strcmp(obj.valueAdjustment,'standard'))
                Qvals = obj.Q(s,:);
            elseif (strcmp(obj.valueAdjustment,'UCB'))
                actionCount = obj.actionCount(s,:);
                actionCount(actionCount==0) = .00001; % for numerical stability
                Qvals = obj.Q(s,:) + obj.c*sqrt(log(obj.totalTime)./actionCount);
            else
                error('obj.valueAdjustment not recognized')
            end
            if (strcmp(obj.actionSelectionMethod,'e-greedy'))
                [maxVal,maxInd] = max(Qvals);
                if (rand < obj.epsilon)
                    action = ceil(rand*length(Qvals));
                else
                    if (length(unique(Qvals)) < length(Qvals)) % for tie breaking
                        maxInds = find(Qvals == max(Qvals));
                        action = maxInds(ceil(rand*length(maxInds)));
                    else
                        action = maxInd;
                    end
                end
            elseif (strcmp(obj.actionSelectionMethod,'softmax'))
                [softVals,action] = mySoftmax(Qvals,obj.tempFunc(obj.totalTime));
            end
            obj.actionCount(s,action) = obj.actionCount(s,action) + 1;
        end
        
        function updateQ(obj,s1,a,r,s2)
            if (strcmp(obj.updateMethod,'q-learning'))
                obj.Q(s1,a) = obj.Q(s1,a) + obj.alpha*(r + obj.gamma*max(obj.Q(s2,:)) - obj.Q(s1,a));
            else
                error('updateQ() is unsure of updateMethod')
            end
        end
        
        function resetActionCount(obj)
            obj.actionCount = zeros(3,2);
        end
     
        function prob = getActionProb(obj,action)
            s = obj.state;
            obj.totalTime = obj.totalTime+1;
            if (strcmp(obj.valueAdjustment,'standard'))
                Qvals = obj.Q(s,:);
            elseif (strcmp(obj.valueAdjustment,'UCB'))
                actionCount = obj.actionCount(s,:);
                actionCount(actionCount==0) = .00001; % for numerical stability
                Qvals = obj.Q(s,:) + obj.c*sqrt(log(obj.totalTime)./actionCount);
            else
                error('obj.valueAdjustment not recognized')
            end
            if (strcmp(obj.actionSelectionMethod,'e-greedy'))
                [maxVal,maxInd] = max(Qvals);
                if (maxInd == action)
                    prob = 1 - obj.epsilon;
                else
                    prob = obj.epsilon;
                end
            elseif (strcmp(obj.actionSelectionMethod,'softmax'))
                [softVals,~] = mySoftmax(Qvals,obj.tempFunc(obj.totalTime));
                prob = softVals(action);
            end
            obj.actionCount(s,action) = obj.actionCount(s,action) + 1;
        end
    end
end
