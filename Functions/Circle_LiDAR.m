function windspeed = Circle_LiDAR(Lidar_x,Lidar_y,Height,num_samples)
% CircleLiDAR parameters:
%     Scanning Pattern: CW, Snapshot
%     Measure Range: 10-550m
%     Weight Function: W
%     Addition: Calculate Helix Center

% Function arguments:
%     Distance x
%     Distance y
%     Height z
%     Discretize num_samples

% FUCK ME! x is pointing at the mother-fucking downwind position

% Derive scanning area
radius = 120;     % assume perfect LiDAR DIEA=240
y0 = Lidar_y;           % center
z0 = Height;            % hub height
measureCenter = [0;0;z0];

% LIDAR sampling --- we want a even circle
[y, z] = meshgrid(linspace(y0 - radius, y0 + radius, num_samples), ...
    linspace(z0 - radius, z0 + radius, num_samples));
mask = (y - y0).^2 + (z - z0).^2 <= radius^2;
y = y(mask);
z = z(mask);

sampled_data = [];
for j = 1:numel(y)
    windspeed_ptr = libpointer('doublePtr', zeros(3,1));
    calllib('QBladeDLL', 'getWindspeed', Lidar_x, y(j), z(j), windspeed_ptr)
    x_vel = windspeed_ptr.value(1); % x component of the windspeed
    y_vel = windspeed_ptr.value(2); % y component of the windspeed
    z_vel = windspeed_ptr.value(3); % z component of the windspeed

    % Projection streamwise speed on LOS to get u_los
    original_los_vector = [Lidar_x; y(j); z(j)] - measureCenter; % los direction vector
    los_unit_vector = original_los_vector/norm(original_los_vector); %unit direction vector
    vel = [x_vel; y_vel; z_vel]; % streamwise speed vector
    u_LOS = dot(vel,los_unit_vector); % should be only one magnitude

    sampled_data(end+1, :) = [Lidar_x, y(j), z(j), x_vel, y_vel, z_vel, u_LOS];
end

% Create structure
windspeed = struct('x', sampled_data(:,1), ...
    'y', sampled_data(:,2), ...
    'z', sampled_data(:,3), ...
    'u_x', sampled_data(:,4), ...
    'u_y', sampled_data(:,5), ...
    'u_z', sampled_data(:,6), ...
    'u_los', sampled_data(:,7));
end