%% Load data
close all
addpath('.\Functions');
windspeed_basic = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_600s_Parallel_Center.mat");
windspeed_turb = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_NTM-A_600s_Parallel_Center.mat");
windspeed_changeTilt = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_600s_Parallel_changeMTilt.mat");
windspeed_changeYaw = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_600s_Parallel_changeMyaw.mat");

dataLiDAR_basic = windspeed_basic.LiDAR_data;
dataLiDAR_turb = windspeed_turb.LiDAR_data;
dataLiDAR_changeTilt = windspeed_changeTilt.LiDAR_data;
dataLiDAR_changeYaw = windspeed_changeYaw.LiDAR_data;

%% Visualization Turbulence
theta = linspace(0, 2*pi, 50);
y_1Dref = 0 + 120 * cos(theta);
z_1Dref = 150 + 120 * sin(theta);

videoFile = ".\Data\Base_vs_Turb.avi";
v = VideoWriter(videoFile);
open(v);

figure('Position', [10, 10, 700, 310]);
for counter = 1:1:data_length(1)  
    subplot(1, 2, 1)
    snapshot = dataLiDAR_basic(counter);
    u_los = snapshot.u_los;
    y = snapshot.y;
    z = snapshot.z;
    wakeCenter = HelixCenter(snapshot, Uin);
    scatter(y, z, 10, u_los, 'filled');
    hold on
    scatter(wakeCenter(1), wakeCenter(2),'red');
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('Base Case', counter)
    colorbar;
    clim([4 8])

    subplot(1, 2, 2)
    snapshot2 = dataLiDAR_turb(counter);
    u_los = snapshot2.u_los;
    y = snapshot2.y;
    z = snapshot2.z;
    wakeCenter2 = HelixCenter(snapshot2, Uin);
    scatter(y, z, 10, u_los, 'filled');
    hold on
    scatter(wakeCenter2(1), wakeCenter2(2),'red');
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('Turbulence', counter)
    colorbar;
    clim([4 8])

    frame = getframe(gcf);
    writeVideo(v, frame);
end 

%% Visualization change tilt
theta = linspace(0, 2*pi, 50);
y_1Dref = 0 + 120 * cos(theta);
z_1Dref = 150 + 120 * sin(theta);

videoFile = ".\Data\Base_vs_Tilt.avi";
v = VideoWriter(videoFile);
open(v);

figure('Position', [10, 10, 700, 310]);
for counter = 1:1:data_length(1)  
    subplot(1, 2, 1)
    snapshot = dataLiDAR_basic(counter);
    u_los = snapshot.u_los;
    y = snapshot.y;
    z = snapshot.z;
    wakeCenter = HelixCenter(snapshot, Uin);
    scatter(y, z, 10, u_los, 'filled');
    hold on
    scatter(wakeCenter(1), wakeCenter(2),'red');
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('Base Case', counter)
    colorbar;
    clim([4 8])

    subplot(1, 2, 2)
    snapshot2 = dataLiDAR_changeTilt(counter);
    u_los = snapshot2.u_los;
    y = snapshot2.y;
    z = snapshot2.z;
    wakeCenter2 = HelixCenter(snapshot2, Uin);
    scatter(y, z, 10, u_los, 'filled');
    hold on
    scatter(wakeCenter2(1), wakeCenter2(2),'red');
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('Change M_{tilt}', counter)
    colorbar;
    clim([4 8])

    frame = getframe(gcf);
    writeVideo(v, frame);
end 

%% Visualization change yaw
theta = linspace(0, 2*pi, 50);
y_1Dref = 0 + 120 * cos(theta);
z_1Dref = 150 + 120 * sin(theta);

videoFile = ".\Data\Base_vs_Yaw.avi";
v = VideoWriter(videoFile);
open(v);

figure('Position', [10, 10, 700, 310]);
for counter = 1:1:data_length(1)  
    subplot(1, 2, 1)
    snapshot = dataLiDAR_basic(counter);
    u_los = snapshot.u_los;
    y = snapshot.y;
    z = snapshot.z;
    wakeCenter = HelixCenter(snapshot, Uin);
    scatter(y, z, 10, u_los, 'filled');
    hold on
    scatter(wakeCenter(1), wakeCenter(2),'red');
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('Base Case', counter)
    colorbar;
    clim([4 8])

    subplot(1, 2, 2)
    snapshot2 = dataLiDAR_changeYaw(counter);
    u_los = snapshot2.u_los;
    y = snapshot2.y;
    z = snapshot2.z;
    wakeCenter2 = HelixCenter(snapshot2, Uin);
    scatter(y, z, 10, u_los, 'filled');
    hold on
    scatter(wakeCenter2(1), wakeCenter2(2),'red');
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('Change M_{yaw}', counter)
    colorbar;
    clim([4 8])

    frame = getframe(gcf);
    writeVideo(v, frame);
end 