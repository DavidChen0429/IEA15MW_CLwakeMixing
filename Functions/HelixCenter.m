function wakeCenter = HelixCenter(snapshot, Uin)
    Original_u_los = snapshot.u_los;
    mean_u_los = mean(Original_u_los);
    threshold = Uin - mean_u_los - 1;                   % threshold
    snapshot.u_losProcessed = Original_u_los - mean_u_los;   % processing
    y = snapshot.y;
    z = snapshot.z;

    % Find indices method1: threshold
    indices = find(snapshot.u_losProcessed >= threshold & sqrt(y.^2 + (z-150).^2) <= 120);
    buf = size(indices);
    if buf(1) == 0
        filtered_indices = sqrt(y.^2 + (z-150).^2) < 120;
        filtered_u_losProcessed = snapshot.u_losProcessed(filtered_indices);
        [~, max_index_within_filtered] = max(filtered_u_losProcessed);
        original_indices = find(filtered_indices);
        indices = original_indices(max_index_within_filtered);
    end

    % Find indices method2: Max with in rotor disc
%     filtered_indices = sqrt(y.^2 + (z-150).^2) < 120;
%     filtered_u_losProcessed = snapshot.u_losProcessed(filtered_indices);
%     [~, max_index_within_filtered] = max(filtered_u_losProcessed);
%     original_indices = find(filtered_indices);
%     indices = original_indices(max_index_within_filtered);

    selected_y = y(indices);
    selected_z = z(indices);
    wakeCenter = [mean(selected_y), mean(selected_z)];
end