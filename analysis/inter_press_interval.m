function [IPI, PRtrials, Incomplete, AbortedPress] = inter_press_interval(SessionData)

% [IPI, PRTRIALS, INCOMPLETE, ABORTEDPRESS] = inter_press_interval(SESSIONDATA)
% 
% This function takes Bpod session data as an input and returns:
%   IPI - A nTrial x 100 matrix of inter-press-intervals for each trial (rows) and press (columns)
%   PRTRIALS - A vector denoting PR trials with a 1
%   INCOMPLETE - A vector denoting incomplete trials with a 1
%   ABORTEDPRESS - A vector denoting if a trial was incomplete and at what press the mouse failed


nTrials = SessionData.nTrials;
IPI = nan(nTrials,100); 
PRtrials = zeros(nTrials,1);
Incomplete = zeros(nTrials,1);
AbortedPress = zeros(nTrials,1);

for i = 1:nTrials
    states = SessionData.RawEvents.Trial{i}.States;
    stateNames = fieldnames(states);
    avoid_overwrite = 0; % for recording AbortedPress
    
    % Find the number of presses required on that trial
    if ~isnan(states.Side_L)
        pressReq = SessionData.TrialSettings(i).GUI.NumPressRequired_L;
        PRtrials(i) = 1;
    elseif ~isnan(states.Side_S)
        pressReq = SessionData.TrialSettings(i).GUI.NumPressRequired_S;
    end
    
    for p = 2:pressReq % skip 1 since there is no press before the first
        structPlace = 3 + p * 2; % position of the relevant struct in 'states'
                                 % add 3 to skip first 3 states:
                                 %     InitialDelay
                                 %     WaitForChoice
                                 %     Side_S/L
                                 % multiply by 2 because each press has 2
                                 % states (wait and debounce)
                                 
        % pull out the time for that press
        leverTimes = getfield(states,[strjoin(stateNames(structPlace))]);
        
        % if the press exists
        if ~isnan(leverTimes(1))
            
            % find the previous lever time and subtract the two
            leverTimes_previous = getfield(states,[strjoin(stateNames(structPlace-2))]);
            IPI(i,p) = leverTimes(end) - leverTimes_previous(end);
            
        % if there are fewer presses than the requirement
        else
            Incomplete(i) = 1;
            
            % do not overwrite a previous lever press
            if avoid_overwrite == 0
                
                % record the last press attempted
                AbortedPress(i) = p;
                avoid_overwrite = 1;
            end
        end
    end   
end

