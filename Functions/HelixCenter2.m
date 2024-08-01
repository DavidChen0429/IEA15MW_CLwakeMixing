function wakeCenter = HelixCenter2(snapshot, Uin)
    Original_u_los = snapshot.u_los;
    mean_u_los = mean(Original_u_los);
    threshold = Uin+0.5;                                         % threshold
    snapshot.u_losProcessed = Original_u_los - mean_u_los;   % processing
    y = snapshot.y;
    z = snapshot.z;

    % Threshold with in rotor disc
    highSpeedIndices = find(snapshot.u_los >= (threshold) & sqrt(y.^2 + (z-150).^2) <= 120);
    buf = size(highSpeedIndices);
    if buf(1) == 0
        filtered_indices = sqrt(y.^2 + (z-150).^2) < 120;
        filtered_u_Processed = snapshot.u_losProcessed(filtered_indices);
        [~, max_index_within_filtered] = max(filtered_u_Processed);
        original_indices = find(filtered_indices);
        validIndices  = original_indices(max_index_within_filtered);
    else
        lowSpeedThreshold = 5; 
        surroundingRadius = 10;
        validIndices = [];
        for i = 1:length(highSpeedIndices)
            idx = highSpeedIndices(i);
            highSpeedY = y(idx);
            highSpeedZ = z(idx);
            surroundingIndices = find(sqrt((y - highSpeedY).^2 + (z - highSpeedZ).^2) <= surroundingRadius);
            surroundingSpeeds = snapshot.u_los(surroundingIndices);
            if min(surroundingSpeeds) <= lowSpeedThreshold
                validIndices = [validIndices; idx];
            end
        end
    end
    
    validY = y(validIndices);
    validZ = z(validIndices);
    wakeCenter = [mean(validY), mean(validZ)];
end