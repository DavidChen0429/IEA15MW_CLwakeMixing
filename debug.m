sampled_data = [];
y0 = 0;
z0 = 0;
radius = 150;
num_samples = 50;

[y, z] = meshgrid(linspace(y0-radius, y0+radius, num_samples), ...
    linspace(z0-radius, z0+radius, num_samples));

mask = (y - y0).^2 + (z - z0).^2 <= radius^2;
y = y(mask);
z = z(mask);
plot(y, z, "*")

if isempty(y) || isempty(z)
    error(['No points within the specified radius. ' ...
    'Please check your radius and num_samples values.'])
end

for j = 1:numel(y)
%     windspeed_ptr = libpointer('doublePtr', zeros(3,1));
%     calllib('QBladeDLL', 'getWindspeed', x, y(j), z(j), windspeed_ptr)
%     x_vel = windspeed_ptr.value(1); % x component of the windspeed
%     sampled_data(end+1, :) = [x, y(j), z(j), x_vel];
    
end