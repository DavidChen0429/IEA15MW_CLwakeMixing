function [r] = referenceGenerator(simTime,Trigger,endTime,type,mag,plotOption)
    % Create reference for CL Helix control
    % Parameter List
    %       simTime: simulation length 
    %       type: type of reference 
    %             step
    %             ramp
    %             ramp&stop
    %             step&step
    %             customize&step
    %             customize&ramp
    %       mag: corresponding helix magnitude

    r = zeros(simTime, 2);    
    % Helix magnitude (from OL simulation)
    if mag == 2
        reference_magnitude = [5.9435 6.0369];
    elseif mag == 3
        reference_magnitude = [8.6257 8.3827];
    end

    if strcmp(type, 'step')
        % Step
        r(Trigger:end, 1) = reference_magnitude(1)*ones(simTime+1-Trigger, 1);   % z_e
        r(Trigger:end, 2) = reference_magnitude(2)*ones(simTime+1-Trigger, 1);   % y_e
    elseif strcmp(type, 'ramp')
        % Ramp
        reference_slope = reference_magnitude/(endTime-Trigger);
        for tt = Trigger:simTime
            r(tt, 1) = reference_slope(1) * (tt - Trigger);   % z_e ramp signal
            r(tt, 2) = reference_slope(2) * (tt - Trigger);   % y_e ramp signal
        end
    elseif strcmp(type, 'ramp&stop')
        % Ramp and Stop
        reference_slope = reference_magnitude/(endTime-Trigger);
        for tt = Trigger:endTime
            r(tt, 1) = reference_slope(1) * (tt - Trigger);   % z_e ramp signal
            r(tt, 2) = reference_slope(2) * (tt - Trigger);   % y_e ramp signal
        end
        r(endTime:end, 1) = reference_magnitude(1)*ones(simTime+1-endTime, 1);
        r(endTime:end, 2) = reference_magnitude(2)*ones(simTime+1-endTime, 1);
    elseif strcmp(type, 'step&step')
        % Complicate step
        steps = cat(2, ...
            0*ones(1, Trigger), 1*ones(1, (simTime-Trigger)/5), ...
            -2*ones(1, (simTime-Trigger)/5), 2*ones(1, (simTime-Trigger)/5), ...
            -1*ones(1, (simTime-Trigger)/5), 0*ones(1, (simTime-Trigger)/5));
        r(:, 1) = steps;
        r(:, 2) = steps;
    elseif strcmp(type, 'zero')
        r = zeros(simTime, 2); 
    elseif strcmp(type, 'customize&step') 
        customMagnitude = [6 9]; % Do NOT exceed 10
        r(Trigger:end, 1) = customMagnitude(1)*ones(simTime+1-Trigger, 1);   % z_e
        r(Trigger:end, 2) = customMagnitude(2)*ones(simTime+1-Trigger, 1);   % y_e
    elseif strcmp(type, 'customize&ramp')
        customMagnitude = [6 9]; % Do NOT exceed 10
        customSlope = customMagnitude/(endTime-Trigger);
        for tt = Trigger:endTime
            r(tt, 1) = customSlope(1) * (tt - Trigger);   % z_e ramp signal
            r(tt, 2) = customSlope(2) * (tt - Trigger);   % y_e ramp signal
        end
        r(endTime:end, 1) = customMagnitude(1)*ones(simTime+1-endTime, 1);
        r(endTime:end, 2) = customMagnitude(2)*ones(simTime+1-endTime, 1);
    end

    if plotOption == 1
        plot(r)
    end
 end