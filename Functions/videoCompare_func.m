function [] = videoCompare_func(data1, data2, Fs, Fc, D, fileName)
dataLiDAR_A = data1.LiDAR_data;
dataLiDAR_B = data2.LiDAR_data;

% wake center traj
wakeCenterZ1 = data1.FF_helixCenter(:, 1);
wakeCenterY1 = data1.FF_helixCenter(:, 2);
wakeCenterZ2 = data2.FF_helixCenter(:, 1);
wakeCenterY2 = data2.FF_helixCenter(:, 2);

% wakeCenterY1_f = lowpassFilter(wakeCenterY1, Fs, Fc);
% wakeCenterZ1_f = lowpassFilter(wakeCenterZ1, Fs, Fc);
% wakeCenterY2_f = lowpassFilter(wakeCenterY2, Fs, Fc);
% wakeCenterZ2_f = lowpassFilter(wakeCenterZ2, Fs, Fc);

data_length = size(dataLiDAR_A);
interval = 10;

% reference rotor disc
theta = linspace(0, 2*pi, 20);
y_1Dref = 0 + D/2 * cos(theta);
z_1Dref = 90 + D/2 * sin(theta);

% ====== Uncomment if save video
% videoFile = fileName; % ".\Data\NTM-B.avi"
% v = VideoWriter(videoFile);
% open(v);

figure('Position', [10, 10, 800, 310]);
for counter = 1:interval:data_length(1)  
    subplot(1, 2, 1)
    snapshot = dataLiDAR_A(counter);
    y = snapshot.y;
    z = snapshot.z;
    scatter(y, z, 10, snapshot.u_los, 'filled');
    hold on
%     scatter(wakeCenterY1(counter), wakeCenterZ1(counter), 'red');
    plot(y_1Dref, z_1Dref, "k:", 'LineWidth',1);
%     plot(wakeCenterY1_f, wakeCenterZ1_f, "r-", 'LineWidth',0.5);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('LiDAR Wind Speed (sec)', round(counter/interval + 1))
    colorbar;
    clim([4 10])

    subplot(1, 2, 2)
    snapshot2 = dataLiDAR_B(counter);
    y = snapshot2.y;
    z = snapshot2.z;
    scatter(y, z, 10, snapshot2.u_los, 'filled');
    hold on
%     scatter(wakeCenterY2(counter), wakeCenterZ2(counter), 'red');
    plot(y_1Dref, z_1Dref, "k:", 'LineWidth',1);
%     plot(wakeCenterY2_f, wakeCenterZ2_f, "r-", 'LineWidth',0.5);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('LiDAR Wind Speed (sec)', round(counter/interval + 1))
    colorbar;
    clim([4 10])
    pause(0.1);

% ====== Uncomment if save video
%     frame = getframe(gcf);
%     writeVideo(v, frame);
end 
end