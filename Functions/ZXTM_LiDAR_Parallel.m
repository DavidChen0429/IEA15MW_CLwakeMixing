function windspeed = ZXTM_LiDAR_Parallel(Lidar_x,Lidar_y,Height,num_samples)
% ZXTM_LIDAR parameters:
%     Scanning Pattern: CW, Ring
%     Scanning Frequency: 50Hz
%     Half-cone-angle: 15 deg
%     Measure Range: 10-550m
%     Weight Function: W

% Function arguments:
%     Distance x
%     Discretize d 

% FUCK ME! x is pointing at the mother-fucking downwind position

% Derive scanning area
%radius = abs(Lidar_x)*tand(15);
radius = 120;     % assume perfect LiDAR DIEA=240
y0 = Lidar_y;           % center
z0 = Height;            % hub height
measureCenter = [0;0;z0];

% Create weight function
%w = 1;  % keep it simple

% LIDAR sampling
theta = linspace(0, 2*pi, num_samples);
theta(end) = [];    % drop the last since 0 == 2pi
y = y0 + radius * cos(theta);
z = z0 + radius * sin(theta);
x = Lidar_x * ones(1, num_samples-1);

posx_ptr = libpointer('doublePtr', x);
posy_ptr = libpointer('doublePtr', y);
posz_ptr = libpointer('doublePtr', z);
velx_ptr = libpointer('doublePtr', zeros(1, num_samples-1));
vely_ptr = libpointer('doublePtr', zeros(1, num_samples-1));
velz_ptr = libpointer('doublePtr', zeros(1, num_samples-1));
calllib('QBladeDLL', 'getWindspeedArray', posx_ptr, posy_ptr, posz_ptr, ...
    velx_ptr,vely_ptr,velz_ptr, num_samples-1)

% Projection to get u_los
original_los_vectors = [240 * ones(num_point, 1), y, z] - repmat(measureCenter, num_point, 1);
los_magnitudes = sqrt(sum(original_los_vectors.^2, 2)); % Magnitude of each LOS vector
los_unit_vectors = original_los_vectors ./ los_magnitudes; % Normalize to get unit vectors

% Calculate vel_los for each point
velocities = [velx_ptr.value, vely_ptr.value, velz_ptr.value]; % Combine velocities into a single matrix
vel_los = sum(velocities .* los_unit_vectors, 2); % Dot product of velocity vectors and LOS unit vectors

% Create structure
windspeed = struct('x', posx_ptr.value, ...
    'y', posy_ptr.value, ...
    'z', posz_ptr.value, ...
    'u_x', velx_ptr.value, ...
    'u_y', vely_ptr.value, ...
    'u_z', velz_ptr.value, ...
    'u_los', vel_los);
end