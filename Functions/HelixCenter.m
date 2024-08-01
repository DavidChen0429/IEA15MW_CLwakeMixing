function wakeCenter = HelixCenter(snapshot, Uin, threshold_diameter)
    Original_u_los = snapshot.u_los;
    mean_u_los = mean(Original_u_los);
    threshold = Uin;                                         % threshold
    snapshot.u_losProcessed = Original_u_los - mean_u_los;   % processing
    y = snapshot.y;
    z = snapshot.z;
    threshold_radii = threshold_diameter/2;

    % Threshold with in rotor disc
    indices = find(snapshot.u_los >= (threshold) & sqrt(y.^2 + (z-90).^2) <= threshold_radii);
    buf = size(indices);
    if buf(1) == 0
        filtered_indices = sqrt(y.^2 + (z-90).^2) < threshold_radii;
        filtered_u_Processed = snapshot.u_losProcessed(filtered_indices);
        [~, max_index_within_filtered] = max(filtered_u_Processed);
        original_indices = find(filtered_indices);
        indices = original_indices(max_index_within_filtered);
    end

    % Max
%     filtered_indices = sqrt(y.^2 + (z-150).^2) < 120;
%     filtered_u_Processed = snapshot.u_losProcessed(filtered_indices);
%     [~, max_index_within_filtered] = max(filtered_u_Processed);
%     original_indices = find(filtered_indices);
%     indices = original_indices(max_index_within_filtered);

    selected_y = y(indices);
    selected_z = z(indices);
    wakeCenter = [mean(selected_z), mean(selected_y)];
end