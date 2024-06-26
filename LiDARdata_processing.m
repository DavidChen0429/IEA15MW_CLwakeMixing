%clc
%clear
close all

%% Load data
windspeed = load('.\Data\MAT\IEA15_Helix_CCW_Str0.3_U8_Uni_600s_1Dd_1Hz_Circle150_windspeedData.mat');
%windspeed = load('.\Data\MAT\IEA15_Helix_CCW_Str0.3_U7_Uni_300s_2Dd_1Hz_windspeedData.mat');
dataLiDAR= windspeed.LiDAR_data;
data_length = size(dataLiDAR);        % length of snapshot
num_snapshot = 100; % one snapshot at one time
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

%% Visualize the Ring
for counter = 1:1:data_length(1)
    snapshot = dataLiDAR(counter);
    u_los = snapshot.u_los;
    y = snapshot.y;
    z = snapshot.z;
    scatter(y, z, 10, u_los, 'filled');
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('LiDAR Wind Speed', counter)
    colorbar;
    clim([3 9])

    filename = sprintf('Data/Figures/LiDAR/1/figure_%d.png', counter);
    saveas(gcf, filename);
end 

%% Average over time (Time Domain)
% get all single point value over simulation time
u_los_average = arrayfun(@(x) x.u_los(1), dataLiDAR);
u_x_average = arrayfun(@(x) x.u_x(1), dataLiDAR);
u_y_average = arrayfun(@(x) x.u_y(1), dataLiDAR);
u_z_average = arrayfun(@(x) x.u_z(1), dataLiDAR);
for counter = 1:1:lengthPoint
    u_los_average = u_los_average + arrayfun(@(x) x.u_los(counter), dataLiDAR);
    u_x_average = u_x_average + arrayfun(@(x) x.u_x(counter), dataLiDAR);
    u_y_average = u_y_average + arrayfun(@(x) x.u_y(counter), dataLiDAR);
    u_z_average = u_z_average + arrayfun(@(x) x.u_z(counter), dataLiDAR);
end
u_x_average = u_x_average / lengthPoint;
u_y_average = u_y_average / lengthPoint;
u_z_average = u_z_average / lengthPoint;
u_los_average = u_los_average / lengthPoint;

figure()
plot(u_los_average)
hold on
plot(u_x_average)
plot(u_y_average)
plot(u_z_average)
% plot(refSine, '--')
xlabel("Time (s)");
ylabel('Speed (m/s)')
legend('u_{los}', 'u_x', 'u_y', 'u_z', 'ref')
title("u_{inflow} Average in Time Domain")
text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
annotation_text = sprintf('Location 2D');
text('Units', 'normalized', 'Position', text_position, ...
    'String', annotation_text, 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'bottom');

%% Average over time (Frequency Domain)
% First, filtered out time for the wake to travel to sampling position
u_los_averageFiltered = u_los_average(100:end);  
Fs = 1;                    % sampling frequency 1Hz
T = 1/Fs;
L = length(u_los_averageFiltered);
t = (0:L-1)*T;  % time vector
n = L;
%n = 2^nextpow2(L);         % Next prime2 number larger than L for faster
f = Fs/n * (0:(n/2));
point_u_los_FFT = fft(u_los_averageFiltered);
P2 = abs(point_u_los_FFT/n);    % normalization
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure()
plot(f,P1,"LineWidth",1);
title("Filtered u_{LOS} Single Point in Frequency Domain")
xlabel("f (Hz)")
ylabel("Magnitude")
text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
annotation_text = sprintf('Location 2D');
text('Units', 'normalized', 'Position', text_position, ...
    'String', annotation_text, 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'bottom');

%% Single point data analysis (Freq, Phase)
% u_x for a single point during simulation intervials
num_point = 8;  % one point over time top point has same phase
pos_x = arrayfun(@(x) x.x(num_point), dataLiDAR);
pos_y = arrayfun(@(x) x.y(num_point), dataLiDAR);
pos_z = arrayfun(@(x) x.z(num_point), dataLiDAR);
point_u_los = arrayfun(@(x) x.u_los(num_point), dataLiDAR);
point_u_x = arrayfun(@(x) x.u_x(num_point), dataLiDAR);
point_u_y = arrayfun(@(x) x.u_y(num_point), dataLiDAR);
point_u_z = arrayfun(@(x) x.u_z(num_point), dataLiDAR);

% Original Signal Plot
figure()
plot(point_u_los)
hold on
plot(point_u_x)
plot(point_u_y)
plot(point_u_z)
plot(refSine, '--')
xlabel("Time (s)");
ylabel('Speed (m/s)')
legend('u_{los}', 'u_x', 'u_y', 'u_z', 'ref')
title("u_{inflow} Single Point in Time Domain")
text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
annotation_text = sprintf('Location 1.86D, point coordinate (%.2f, %.2f)', pos_y(1), pos_z(1));
text('Units', 'normalized', 'Position', text_position, ...
    'String', annotation_text, 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'bottom');

%%
% plot(arrayfun(@(x) x.u_los(1), dataLiDAR))
% hold on
% plot(arrayfun(@(x) x.u_los(3), dataLiDAR))
% legend('1','3')

for i = 1:10
    point_u_los = arrayfun(@(x) x.u_los(i), dataLiDAR);
    plot(point_u_los)
    hold on
end

%% FFT Signal point overtime 
point_u_losFiltered = point_u_los(100:end);   % filter out nonexcited period
Fs = 1;                    % sampling frequency 1Hz
T = 1/Fs;
L = length(point_u_losFiltered);
n = 2^nextpow2(L);         % Next prime2 number larger than L for faster
f = Fs/n * (0:(n/2));
point_u_los_FFT = fft(point_u_losFiltered, n);
P2 = abs(point_u_los_FFT/n);    % normalization
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure()
plot(f,P1,"LineWidth",1);
title("Filtered u_{LOS} Single Point in Frequency Domain")
xlabel("f (Hz)")
ylabel("Magnitude")
text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
annotation_text = sprintf('Location 2D');
text('Units', 'normalized', 'Position', text_position, ...
    'String', annotation_text, 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'bottom');

%% Turbulence Analysis
% FFT whole ring in one snapshot (Turbulence) 
% Filtered from 50
for i = 1:1:data_length(1)
    u_los = dataLiDAR(i).u_los;
    L2 = length(u_los);
    snapshot_FFT = fft(u_los);
    scatter(dataLiDAR(i).y, dataLiDAR(i).z, abs(snapshot_FFT), 'filled');
    xlabel('Y [m]')
    ylabel('Z [m]')
    title('LiDAR FFT for snapshot', i)
    colorbar;
    filename = sprintf('Data/Figures/Version_ringFFT/figureFFT_%d.png', i);
    saveas(gcf, filename);
end