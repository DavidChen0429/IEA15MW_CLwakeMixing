function [] = videoCompare_func(data1, data2, Fs, Fc)
dataLiDAR_A = data1.LiDAR_data;
dataLiDAR_B = data2.LiDAR_data;

% wake center traj
wakeCenterY1 = arrayfun(@(x) x.centerY, dataLiDAR_A);
wakeCenterZ1 = arrayfun(@(x) x.centerZ, dataLiDAR_A);
wakeCenterY2 = arrayfun(@(x) x.centerY, dataLiDAR_B);
wakeCenterZ2 = arrayfun(@(x) x.centerZ, dataLiDAR_B);

wakeCenterY1_f = lowpassFilter(wakeCenterY1, Fs, Fc);
wakeCenterZ1_f = lowpassFilter(wakeCenterZ1, Fs, Fc);
wakeCenterY2_f = lowpassFilter(wakeCenterY2, Fs, Fc);
wakeCenterZ2_f = lowpassFilter(wakeCenterZ2, Fs, Fc);

data_length = size(dataLiDAR_A);
interval = 10;

% reference rotor disc
theta = linspace(0, 2*pi, 20);
y_1Dref = 0 + 120 * cos(theta);
z_1Dref = 150 + 120 * sin(theta);

% ====== Uncomment if save video
% videoFile = ".\Data\Me_tilt_l.avi";
% v = VideoWriter(videoFile);
% open(v);

figure('Position', [10, 10, 700, 310]);
for counter = 1:interval:data_length(1)  
    subplot(1, 2, 1)
    snapshot = dataLiDAR_A(counter);
    y = snapshot.y;
    z = snapshot.z;
    scatter(y, z, 10, snapshot.u_los, 'filled');
    hold on
    scatter(snapshot.centerY, snapshot.centerZ,'red');
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
    plot(wakeCenterY1_f, wakeCenterZ1_f, "r-", 'LineWidth',0.5);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('LiDAR Wind Speed (sec)', round(counter/interval + 1))
    colorbar;
    clim([4 8])

    subplot(1, 2, 2)
    snapshot2 = dataLiDAR_B(counter);
    y = snapshot2.y;
    z = snapshot2.z;
    scatter(y, z, 10, snapshot2.u_los, 'filled');
    hold on
    scatter(snapshot2.centerY, snapshot2.centerZ,'red');
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
    plot(wakeCenterY2_f, wakeCenterZ2_f, "r-", 'LineWidth',0.5);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('LiDAR Wind Speed (sec)', round(counter/interval + 1))
    colorbar;
    clim([4 8])
    pause(0.1);

%     frame = getframe(gcf);
%     writeVideo(v, frame);
end 
end

