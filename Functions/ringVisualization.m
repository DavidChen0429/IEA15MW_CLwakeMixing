function [] = ringVisualization(dataLiDAR, D)
interval = 10;
data_length = size(dataLiDAR);
% Rotor Disc
theta = linspace(0, 2*pi, 20);
y_1Dref = 0 + D/2 * cos(theta);
z_1Dref = 90 + D/2 * sin(theta);

% Visualization
figure;
for counter = 1:interval:data_length(1)  
    snapshot = dataLiDAR(counter);
    u_los = snapshot.u_los;
    y = snapshot.y;
    z = snapshot.z;
    scatter(y, z, 10, u_los, 'filled');

%     u_x = snapshot.u_x;
%     u_y = snapshot.u_y;
%     u_z = snapshot.u_z;
%     magnitude_speed = sqrt(u_x.^2 + u_y.^2 + u_z.^2);
%     scatter(y, z, 10, magnitude_speed, 'filled');

    hold on
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
%     scatter(snapshot.centerY, snapshot.centerZ,'red');
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('LiDAR Wind Speed (sec)', round(counter/interval + 1))
    colorbar;
    clim([4 11])
    pause(0.1);
end
end