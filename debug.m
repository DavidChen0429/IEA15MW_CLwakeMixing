y0 = 0;
z0 = 150;
radius = 120;
num_samples = 100;
measureCenter = [0 0 z0];

[y, z] = meshgrid(linspace(y0-radius, y0+radius, num_samples), ...
    linspace(z0-radius, z0+radius, num_samples));

mask = (y - y0).^2 + (z - z0).^2 <= radius^2;
y = y(mask);
z = z(mask);
plot(y, z, "*")
buf = size(y);
num_point = buf(1)

if isempty(y) || isempty(z)
    error(['No points within the specified radius. ' ...
    'Please check your radius and num_samples values.'])
end

original_los_vectors = [240 * ones(num_point, 1), y, z] - repmat(measureCenter, num_point, 1);
los_magnitudes = sqrt(sum(original_los_vectors.^2, 2)); % Magnitude of each LOS vector
los_unit_vectors = original_los_vectors ./ los_magnitudes; % Normalize to get unit vectors

%%
radius = 120;
y0 = 0;
z0 = 150;
theta = linspace(0, 2*pi, 100);
y = y0 + radius * cos(theta);
z = z0 + radius * sin(theta);
plot(y, z, "-", 'LineWidth',3);