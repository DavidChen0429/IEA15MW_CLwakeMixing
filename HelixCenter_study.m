%clc
%clear
close all
addpath('.\Functions');

%% Load data
windspeed = load('.\Data\MAT\LiDAR_sampling\IEA15_Helix_CCW_Str0.3_U8_Uni_600s_1Dd_1Hz_Circle150_windspeedData');
% windspeed = load('.\Data\MAT\LiDAR_sampling\IEA15_Helix_CCW_Str0.3_U8_Uni_300s_1Dd_1Hz_Circle_windspeedData.mat');

dataLiDAR= windspeed.LiDAR_data;
data_length = size(dataLiDAR);        % length of snapshot
lengthPoint = length(dataLiDAR(1).x);
measurementPos = 450;
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

%% Helix Center Study
theta = linspace(0, 2*pi, 20);
y_1Dref = 0 + 120 * cos(theta);
z_1Dref = 150 + 120 * sin(theta);

snapshot = dataLiDAR(299);
y = snapshot.y;
z = snapshot.z;
u_losOriginal = snapshot.u_los;
u_losProcessed = u_losOriginal - mean(u_losOriginal);
wakeCenter = HelixCenter(snapshot, Uin);

% Visualization
scatter(y, z, 10, u_losProcessed, 'filled');
hold on
plot(y_1Dref, z_1Dref, "k-", 'LineWidth',2);
scatter(wakeCenter(1), wakeCenter(2),'red');
xlabel('Y [m]')
ylabel('Z [m]')
title('LiDAR Wind Speed')
colorbar;
upperbound = 8 - mean(u_losOriginal);
lowerbound = 4 - mean(u_losOriginal);
clim([lowerbound upperbound])

%% Calculate the Helix center
wake_center = [];
time = 300;
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
save('300sMix1874.mat', 'wake_center');


%% Visualize the Helix Center and Helix
% Reference 1D ring
theta = linspace(0, 2*pi, 50);
y_1Dref = 0 + 120 * cos(theta);
z_1Dref = 150 + 120 * sin(theta);
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
    pause(0.1);
end 

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

%% Compare Different Resolution 
a = load(".\Data\MAT\Helix_wake_center\300sMix1874.mat");
b = load(".\Data\MAT\Helix_wake_center\300sMix276.mat");
figure();
sgtitle('Different Resolution');
subplot(2,1,1);
plot(t,a.wake_center(:,1))
hold on 
plot(t,b.wake_center(:,1))
ylim([-120 120]);
xlabel('Time [s]')
ylabel('Y [m]')
legend('1874', '276')
subplot(2,1,2);
plot(t,a.wake_center(:,2))
hold on 
plot(t,b.wake_center(:,2))
ylim([30 270]);
xlabel('Time [s]')
ylabel('Z [m]')
legend('1874', '276')