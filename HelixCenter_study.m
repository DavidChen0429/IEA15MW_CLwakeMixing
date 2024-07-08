%clc
%clear
close all
addpath('.\Functions');

%% Load data
% windspeed = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_600s_Parallel_changeMTilt.mat");
windspeed = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_600s_Parallel_Center.mat");

dataLiDAR= windspeed.LiDAR_data;
data_length = size(dataLiDAR);        % length of snapshot
lengthPoint = length(dataLiDAR(1).x);
measurementPos = 240;
Uin = 8;
Str = 0.3;           % Strouhal number 
DIEA15 = 240;
Freq = Str*Uin/DIEA15;      % From Str, in Hz

% Reference signal
timeWakeTravel0 = round(measurementPos/Uin);
timeWakeTravel0 = 45;
t = linspace(0, data_length(1)-timeWakeTravel0, data_length(1)-timeWakeTravel0);
refSine1 = 2 * sin(2*pi*Freq*t);
bufferSine = zeros(1, timeWakeTravel0);
refSine = [bufferSine, refSine1] + 6.5;

% snapshot variables: 
%   Position info:          x,y,z
%   Streamwise speed info:  u_x,u_y,u_z
%   LOS speed info:         u_los

%% Pure visualization
theta = linspace(0, 2*pi, 20);
y_1Dref = 0 + 120 * cos(theta);
z_1Dref = 150 + 120 * sin(theta);

% Visualization
figure;
for counter = 1:1:data_length(1)  
    snapshot = dataLiDAR(counter);
    u_los = snapshot.u_los;
    y = snapshot.y;
    z = snapshot.z;
%     scatter(y, z, 10, u_los, 'filled');

    u_x = snapshot.u_x;
    u_y = snapshot.u_y;
    u_z = snapshot.u_z;
    magnitude_speed = sqrt(u_x.^2 + u_y.^2 + u_z.^2);
    scatter(y, z, 10, magnitude_speed, 'filled');

    hold on
    plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
    hold off;
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('LiDAR Wind Speed', counter)
    colorbar;
    clim([4 9])
    pause(0.1);
end

%% Calculate the Helix center
wake_center = [];
time = 600;
t = linspace(1, time, time);
for counter = 1:1:time
    snapshot = dataLiDAR(counter);
    wakeCenter = HelixCenter(snapshot, Uin);
    wake_center(end+1, :) = [wakeCenter(1), wakeCenter(2)];
end 
figure();
sgtitle('Helix Wake Center');
subplot(2,1,1);
plot(t, wake_center(:,1), 'r-');
ylim([-120 120]);
xlabel('Time [s]')
ylabel('Y [m]')
subplot(2,1,2);
plot(t, wake_center(:,2), 'b-');
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')
save('.\Data\MAT\Helix_wake_center\600sMix2724_changeTilt.mat', 'wake_center');

%% Compare the real-time computation
% Offline computation
wake_center = [];
time = 600;
t = linspace(1, time, time);
for counter = 1:1:time
    snapshot = dataLiDAR(counter);
    wakeCenter = HelixCenter(snapshot, Uin);
    wake_center(end+1, :) = [wakeCenter(1), wakeCenter(2)];
end 

% Online computation
wakeCenterY = arrayfun(@(x) x.centerY, dataLiDAR);
wakeCenterZ = arrayfun(@(x) x.centerZ, dataLiDAR);

figure();
sgtitle('Helix Wake Center');
subplot(2,1,1);
plot(t, wake_center(:,1), 'r-');
hold on
plot(t, wakeCenterY)
legend('offline', 'online')
ylim([-120 120]);
xlabel('Time [s]')
ylabel('Y [m]')
subplot(2,1,2);
plot(t, wake_center(:,2), 'b-');
hold on
plot(t, wakeCenterZ)
legend('offline', 'online')
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')

%% Visualize the Helix Center and Helix
% Reference 1D ring
theta = linspace(0, 2*pi, 50);
y_1Dref = 0 + 120 * cos(theta);
z_1Dref = 150 + 120 * sin(theta);

figure('Position', [10, 10, 700, 310]);
for counter = 1:1:data_length(1)  
    subplot(1, 2, 1)
    snapshot = dataLiDAR(counter);
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
    title('LiDAR Wind Speed', counter)
    colorbar;
    clim([4 8])

    subplot(1, 2, 2)
    snapshot = dataLiDAR(counter);
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
    title('LiDAR Wind Speed', counter)
    colorbar;
    clim([4 8])
    pause(0.1);
end 


%% Save as video
% Reference 1D ring
theta = linspace(0, 2*pi, 50);
y_1Dref = 0 + 120 * cos(theta);
z_1Dref = 150 + 120 * sin(theta);

videoFile = ".\Data\Helix_basic.avi";
v = VideoWriter(videoFile);
open(v);

figure;
for counter = 1:1:data_length(1)
    snapshot = dataLiDAR(counter);
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
    title('LiDAR Wind Speed', counter)
    colorbar;
    clim([4 8])
%     pause(0.1);

    frame = getframe(gcf);
    writeVideo(v, frame);
end 

close(v);

%% Compare two methods
a = load(".\Data\MAT\Helix_wake_center\300sThreshold1874.mat");
b = load(".\Data\MAT\Helix_wake_center\300sMax1874.mat");
figure();
sgtitle('Different Method of Helix Center');
subplot(2,1,1);
plot(t,a.wake_center(:,1))
hold on 
plot(t,b.wake_center(:,1))
ylim([-120 120]);
xlabel('Time [s]')
ylabel('Y [m]')
legend('threshold', 'max')
subplot(2,1,2);
plot(t,a.wake_center(:,2))
hold on 
plot(t,b.wake_center(:,2))
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')
legend('threshold', 'max')

%% Compare u_los and u_mag
a = load(".\Data\MAT\Helix_wake_center\Turb600sMix2724_magU.mat");
b = load(".\Data\MAT\Helix_wake_center\Turb600sMix2724_ulos.mat");
c = load(".\Data\MAT\Helix_wake_center\600sMix2724_magU.mat");
d = load(".\Data\MAT\Helix_wake_center\600sMix2724_ulos.mat");
figure();
sgtitle('Different u of computing Helix Center');
subplot(2,1,1);
plot(t,a.wake_center(:,1))
hold on 
plot(t,b.wake_center(:,1))
plot(t,c.wake_center(:,1))
plot(t,d.wake_center(:,1))
ylim([-120 120]);
xlabel('Time [s]')
ylabel('Y [m]')
legend('|U|', 'U_{los}','|U| turb', 'U_{los} turb')
subplot(2,1,2);
plot(t,a.wake_center(:,2))
hold on 
plot(t,b.wake_center(:,2))
plot(t,c.wake_center(:,2))
plot(t,d.wake_center(:,2))
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')
legend('|U|', 'U_{los}','|U| turb', 'U_{los} turb')

%% Compare Different Resolution 
% Uniform case
a = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point1874_Timestep0.1_600s_Parallel_Center.mat");
b = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point1184_Timestep0.1_600s_Parallel_Center.mat");
c = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point648_Timestep0.1_600s_Parallel_Center.mat");
d = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point276_Timestep0.1_600s_Parallel_Center.mat");
e = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_600s_Parallel_Center.mat");
f = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point60_Timestep0.1_600s_Parallel_Center.mat");
a_wakeCenterY = arrayfun(@(x) x.centerY, a.LiDAR_data);
a_wakeCenterZ = arrayfun(@(x) x.centerZ, a.LiDAR_data);
b_wakeCenterY = arrayfun(@(x) x.centerY, b.LiDAR_data);
b_wakeCenterZ = arrayfun(@(x) x.centerZ, b.LiDAR_data);
c_wakeCenterY = arrayfun(@(x) x.centerY, c.LiDAR_data);
c_wakeCenterZ = arrayfun(@(x) x.centerZ, c.LiDAR_data);
d_wakeCenterY = arrayfun(@(x) x.centerY, d.LiDAR_data);
d_wakeCenterZ = arrayfun(@(x) x.centerZ, d.LiDAR_data);
e_wakeCenterY = arrayfun(@(x) x.centerY, e.LiDAR_data);
e_wakeCenterZ = arrayfun(@(x) x.centerZ, e.LiDAR_data);
f_wakeCenterY = arrayfun(@(x) x.centerY, f.LiDAR_data);
f_wakeCenterZ = arrayfun(@(x) x.centerZ, f.LiDAR_data);
t = linspace(1, 600, 600);

figure();
sgtitle('Different Resolution Uniform Wind');
subplot(2,1,1);
plot(t, e_wakeCenterY)
hold on 
plot(t, a_wakeCenterY)
plot(t, b_wakeCenterY)
plot(t, c_wakeCenterY)
plot(t, d_wakeCenterY)
plot(t, f_wakeCenterY)
ylim([-120 120]);
xlabel('Time [s]')
ylabel('Y [m]')
legend('2724','1874','1184','648','276','60')
subplot(2,1,2);
plot(t, e_wakeCenterZ)
hold on 
plot(t, a_wakeCenterZ)
plot(t, b_wakeCenterZ)
plot(t, c_wakeCenterZ)
plot(t, d_wakeCenterZ)
plot(t, f_wakeCenterZ)
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')
legend('2724','1874','1184','648','276','60')

% Turbulence case
a2 = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point1874_Timestep0.1_NTM-C_600s_Parallel_Center.mat");
b2 = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point1184_Timestep0.1_NTM-C_600s_Parallel_Center.mat");
c2 = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point648_Timestep0.1_NTM-C_600s_Parallel_Center.mat");
d2 = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point276_Timestep0.1_NTM-C_600s_Parallel_Center.mat");
e2 = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point149_Timestep0.1_NTM-C_600s_Parallel_Center.mat");
f2 = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_NTM-C_600s_Parallel_Center.mat");
a_wakeCenterY2 = arrayfun(@(x) x.centerY, a2.LiDAR_data);
a_wakeCenterZ2 = arrayfun(@(x) x.centerZ, a2.LiDAR_data);
b_wakeCenterY2 = arrayfun(@(x) x.centerY, b2.LiDAR_data);
b_wakeCenterZ2 = arrayfun(@(x) x.centerZ, b2.LiDAR_data);
c_wakeCenterY2 = arrayfun(@(x) x.centerY, c2.LiDAR_data);
c_wakeCenterZ2 = arrayfun(@(x) x.centerZ, c2.LiDAR_data);
d_wakeCenterY2 = arrayfun(@(x) x.centerY, d2.LiDAR_data);
d_wakeCenterZ2 = arrayfun(@(x) x.centerZ, d2.LiDAR_data);
e_wakeCenterY2 = arrayfun(@(x) x.centerY, e2.LiDAR_data);
e_wakeCenterZ2 = arrayfun(@(x) x.centerZ, e2.LiDAR_data);
f_wakeCenterY2 = arrayfun(@(x) x.centerY, f2.LiDAR_data);
f_wakeCenterZ2 = arrayfun(@(x) x.centerZ, f2.LiDAR_data);
t = linspace(1, 600, 600);

figure();
sgtitle('Different Resolution Turbulent Wind');
subplot(2,1,1);
plot(t, f_wakeCenterY2)
hold on 
plot(t, b_wakeCenterY2)
plot(t, c_wakeCenterY2)
plot(t, d_wakeCenterY2)
plot(t, e_wakeCenterY2)
plot(t, a_wakeCenterY2)
ylim([-120 120]);
xlabel('Time [s]')
ylabel('Y [m]')
legend('2724','1874','1184','648','276','149')
subplot(2,1,2);
plot(t, f_wakeCenterZ2)
hold on 
plot(t, b_wakeCenterZ2)
plot(t, c_wakeCenterZ2)
plot(t, d_wakeCenterZ2)
plot(t, e_wakeCenterZ2)
plot(t, a_wakeCenterZ2)
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')
legend('2724','1874','1184','648','276','149')

%% Compare different cases
Tur_NTM_A = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_NTM-A_600s_Parallel_Center.mat");
Tur_NTM_B = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_NTM-B_600s_Parallel_Center.mat");
Tur_NTM_C = load(".\Data\MAT\LiDAR_sampling\Turbulence\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_NTM-C_600s_Parallel_Center.mat");
basecase = load(".\Data\MAT\LiDAR_sampling\Uni\Str0.3_U8_1Dd_1Hz\Point2724_Timestep0.1_600s_Parallel_Center.mat");

TurbA_Y = arrayfun(@(x) x.centerY, Tur_NTM_A.LiDAR_data);
TurbA_Z = arrayfun(@(x) x.centerZ, Tur_NTM_A.LiDAR_data);
TurbB_Y = arrayfun(@(x) x.centerY, Tur_NTM_B.LiDAR_data);
TurbB_Z = arrayfun(@(x) x.centerZ, Tur_NTM_B.LiDAR_data);
TurbC_Y = arrayfun(@(x) x.centerY, Tur_NTM_C.LiDAR_data);
TurbC_Z = arrayfun(@(x) x.centerZ, Tur_NTM_C.LiDAR_data);
base_Y = arrayfun(@(x) x.centerY, basecase.LiDAR_data);
base_Z = arrayfun(@(x) x.centerZ, basecase.LiDAR_data);
t = linspace(1, 600, 600);

figure();
sgtitle('Different Cases');
subplot(2,1,1);
plot(t, TurbA_Y)
hold on 
plot(t, TurbB_Y)
plot(t, TurbC_Y)
plot(t, base_Y)
% ylim([-120 120]);
ylim([-50 50]);
xlabel('Time [s]')
ylabel('Y [m]')
legend('NTM-A','NTM-B','NTM-C','Base')
subplot(2,1,2);
plot(t, TurbA_Z)
hold on 
plot(t, TurbB_Z)
plot(t, TurbC_Z)
plot(t, base_Z)
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')
legend('NTM-A','NTM-B','NTM-C','Base')

%% Compare Mtilt
% Load wake center data
a = load(".\Data\MAT\Helix_wake_center\600sMix2724_basic.mat");
b = load(".\Data\MAT\Helix_wake_center\600sMix2724_changeTilt.mat");

% Mtilt and Myaw
midPoint = 600 / 2;
Helix_amplitude = 4;   
amplitudeTilt = Helix_amplitude * (t <= midPoint) + 2 * Helix_amplitude * (t > midPoint);
sigTilt = amplitudeTilt .* sin(2*pi*Freq*t);  
sigTiltoriginal = Helix_amplitude * sin(2*pi*Freq*t);          
sigYaworiginal = Helix_amplitude * sin(2*pi*Freq*t + pi/2);  % CCW
sigYaw = amplitudeTilt .* sin(2*pi*Freq*t + pi/2);  % CCW

figure();
sgtitle('M_{tilt} and Helix center');
subplot(3,1,1);
plot(t,a.wake_center(:,1),'m')
hold on 
plot(t,b.wake_center(:,1),'b')
hold off
ylim([-120 120]);
xlabel('Time [s]')
ylabel('Y [m]')
legend('original', 'M^{change}_{tilt}')
subplot(3,1,2);
plot(t,a.wake_center(:,2),'m')
hold on 
plot(t,b.wake_center(:,2),'b')
hold off
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')
legend('original', 'M^{change}_{tilt}')
subplot(3,1,3);
plot(t, sigTilt,'b')
hold on
plot(t, sigYaworiginal,'m')
plot(t, sigTiltoriginal, 'm--')
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('M^{change}_{tilt}','M_{yaw}','M_{tilt}')

%% Compare Myaw
% Load wake center data
a = load(".\Data\MAT\Helix_wake_center\600sMix2724_basic.mat");
% b = load(".\Data\MAT\Helix_wake_center\600sMix2724_changTilt.mat");
c = load(".\Data\MAT\Helix_wake_center\600sMix2724_changeYaw.mat");

% Mtilt and Myaw
midPoint = 600 / 2;
Helix_amplitude = 4;   
amplitudeTilt = Helix_amplitude * (t <= midPoint) + 2 * Helix_amplitude * (t > midPoint);
sigTilt = amplitudeTilt .* sin(2*pi*Freq*t);  
sigTiltoriginal = Helix_amplitude * sin(2*pi*Freq*t);          
sigYaworiginal = Helix_amplitude * sin(2*pi*Freq*t + pi/2);  % CCW
sigYaw = amplitudeTilt .* sin(2*pi*Freq*t + pi/2);  % CCW

figure();
sgtitle('M_{yaw} and Helix center');
subplot(3,1,1);
plot(t,a.wake_center(:,1),'m')
hold on 
plot(t,c.wake_center(:,1),'b')
hold off
ylim([-120 120]);
xlabel('Time [s]')
ylabel('Y [m]')
legend('original', 'M^{change}_{yaw}')
subplot(3,1,2);
plot(t,a.wake_center(:,2),'m')
hold on 
plot(t,c.wake_center(:,2),'b')
hold off
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')
legend('original', 'M^{change}_{yaw}')
subplot(3,1,3);
plot(t, sigTiltoriginal,'m')
hold on
plot(t, sigYaworiginal,'m--')
plot(t, sigYaw, 'b')
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('M_{tilt}','M_{yaw}','M^{change}_{yaw}')