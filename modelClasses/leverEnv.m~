classdef leverEnv < matlab.mixin.Copyable
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        itiCost % cost associated with waiting through the inter-trial interval
        leverPressCost % cost associated with pressing a lever
        Sreward % small reward size
        Lreward % large reward size
        agent % associated agent (animal)
        nPressesForReward_S % # of presses at SR required for reward
        nPressesForReward_L % # of presses at LR required for reward (changes)
        nPressesS % current # of consecutive presses at SR
        nPressesL % current # of consecutive presses at LR
        rewardFunc
    end
    
    methods
        function obj = leverEnv(opts,agent)
            % Initialize the environment with input parameters and
            % assocated agent
            obj.itiCost = opts.itiCost;
            obj.leverPressCost = opts.leverPressCost;
            obj.Sreward = opts.Sreward;
            obj.Lreward = opts.Lreward;
            obj.nPressesForReward_S = opts.nPressesForReward_S;
            obj.nPressesForReward_L = opts.nPressesForReward_L;
            obj.agent = agent;
            obj.nPressesS = 0;
            obj.nPressesL = 0;
            obj.rewardFunc = opts.rewardFunc;
        end
        
        function [r,newS] = processAction(obj,action)
            % Process the agent's action and emit a reward and new state
            s = obj.agent.state; % agent's current state
            if (s == 1) % If the current state is a trial start
                newS = s + action;
                r = obj.leverPressCost; % no cost in 'deciding' which lever to choose
            elseif (s == 2) % If the agent's last choice was to push the SR lever
                if (action == 1) % If the current choice is to push the lever
                    obj.nPressesS = obj.nPressesS + 1;
                    if (obj.nPressesS == obj.nPressesForReward_S) % If the SR lever has been pushed enough for reward
                        newS = 1; % Move the agent back to the trial start
                        r = obj.Sreward - obj.leverPressCost - obj.itiCost; % Reward = R(SR) - (cost of lever pressing + cost of waiting for next trial)
                        obj.nPressesS = 0; % reset the # of presses at SR
                    else % If additional SR presses are needed for reward
                        newS = 2; % the state stays the same
                        r = obj.leverPressCost; % the agent incurs the lever press cost
                    end
                elseif (action == 2) % If the agent decideds to abort the current trial
                    newS = 1; % Reset the agent to the trial start
                    r = obj.itiCost; % The agent incurs the iti cost
                    obj.nPressesS = 0; % reset the # of presses at SR
                end
            elseif (s == 3) % If the agent's last choice was to push the LR lever
                if (action == 1) % If the current choice is to push the lever
                    obj.nPressesL = obj.nPressesL + 1;
                    if (obj.nPressesL == obj.nPressesForReward_L) % If the LR lever has been pushed enough for reward
                        newS = 1; % Move the agent back to the trial start
                        r = obj.Lreward - obj.leverPressCost - obj.itiCost; % Reward = R(LR) - (cost of lever pressing + cost of waiting for next trial)
                        obj.nPressesL = 0; % Reset the # of consecutive presses at LR to 0
                        obj.nPressesForReward_L = obj.nPressesForReward_L+1; % Increment the # of presses required for reward at LR by 1
                    else % If more LR presses are required for reward
                        newS = 3; % the state stays the same
                        r = obj.leverPressCost; % the agent incurs the lever press cost
                    end
                elseif (action == 2) % If the agent decides to abort the current trial
                    newS = 1; % Reset the agent to the trial start
                    r = obj.itiCost; % The agent incurs the iti cost
                    obj.nPressesL = 0; % Reset the # of consecutive presses at LR to 0
                end
            end
        end
        
    end
end