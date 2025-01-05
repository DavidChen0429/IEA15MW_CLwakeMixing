function [] = videoCompare_func(data1, data2, D, fileName)
dataLiDAR_A = data1.LiDAR_data;
dataLiDAR_B = data2.LiDAR_data;

data_length = size(dataLiDAR_A);
interval = 10;

% reference rotor disc
theta = linspace(0, 2*pi, 20);
y_1Dref = 0 + D/2 * cos(theta);
z_1Dref = 90 + D/2 * sin(theta);

% ====== Uncomment if save video
videoFile = fileName; % ".\Data\NTM-B.avi"
v = VideoWriter(videoFile);
open(v);

figure('Position', [10, 10, 800, 310]);
for counter = 1:interval:data_length(1)  
    subplot(1, 2, 1)
    snapshot = dataLiDAR_A(counter);
    y = snapshot.y;
    z = snapshot.z;
    scatter(y, z, 10, snapshot.u_los, 'filled');
    hold on
    plot(y_1Dref, z_1Dref, "k:", 'LineWidth',1);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('Baseline (sec)', round(counter/interval + 1))
    colorbarHandle = colorbar;
    ylabel(colorbarHandle, 'u [m/s]');
    clim([4 10])

    subplot(1, 2, 2)
    snapshot2 = dataLiDAR_B(counter);
    y = snapshot2.y;
    z = snapshot2.z;
    scatter(y, z, 10, snapshot2.u_los, 'filled');
    hold on
    plot(y_1Dref, z_1Dref, "k:", 'LineWidth',1);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('CCW Helix (sec)', round(counter/interval + 1))
    colorbarHandle = colorbar;
    ylabel(colorbarHandle, 'u [m/s]');
    clim([4 10])
    pause(0.1);

% ====== Uncomment if save video
    frame = getframe(gcf);
    writeVideo(v, frame);
end 
end