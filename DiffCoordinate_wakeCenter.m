%% This file is created to learn the different cooridnate for the helix center
clear
close all
addpath('.\Functions');

%% Basic information definition
fileName = 'Point2724_300s_Center_helix.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';
windspeed = load([dataPath caseName fileName]);
dataLiDAR= windspeed.LiDAR_data;
filteredIndex = 60;
dataLiDAR = dataLiDAR(filteredIndex:end);
data_length = size(dataLiDAR);        % length of snapshot
lengthPoint = length(dataLiDAR(1).x);

U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;     % Rotor diameter
simTime = 10000;     % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds
Str = 0.3;                          % Strouhal number
Helix_amplitude = 4;                % Helix amplitude                
Freq = Str*U_inflow/D_IEA15MW;      % From Str, in Hz
omega_e = Freq*2*pi;

% Tranform matrix to Helix coordinate frame
% R_helix = [1 0 0; 
%            0 cos(Freq) -sin(Freq); 
%            0 sin(Freq) cos(Freq)];
% inverseR_helix = [1 0 0; 
%                   0 cos(Freq) sin(Freq); 
%                   0 -sin(Freq) cos(Freq)];

%% From Fixed Frame --> Helix Frame, CW
% Fixed Frame
wakeCenterY = arrayfun(@(x) x.centerY, dataLiDAR);
wakeCenterZ = arrayfun(@(x) x.centerZ, dataLiDAR);
t = linspace(1, data_length(1), data_length(1));
% FFT_func(wakeCenterY, 1, 1)

% Store wake center information in helix frame
wakecenterY_helixFrame_store = zeros(data_length(1), 1);
wakecenterZ_helixFrame_store = zeros(data_length(1), 1);
wakecenterY_fixFrame_store = zeros(data_length(1), 1);
wakecenterZ_fixFrame_store = zeros(data_length(1), 1);

for i = 1:1:data_length(1)
    R_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
               -sin(omega_e*t(i)) cos(omega_e*t(i))];
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    centerY = wakeCenterY(i);
    centerZ = wakeCenterZ(i);

    % Fixed Frame ---> Helix Frame
    center_helixFrame = R_helix * [centerZ; centerY]; % tilt, yaw
    wakecenterY_helixFrame_store(i) = center_helixFrame(2);
    wakecenterZ_helixFrame_store(i) = center_helixFrame(1);

    % Helix Frame ---> Fixed Frame
    center_fixedFrame = invR_helix * center_helixFrame; 
    wakecenterY_fixFrame_store(i) = center_fixedFrame(2);
    wakecenterZ_fixFrame_store(i) = center_fixedFrame(1);
end

figure();
subplot(3,1,1);
plot(t, wakeCenterY + mean(wakeCenterZ))
hold on;
plot(t, wakeCenterZ)
hold off;
xlabel('Time [s]')
ylabel('Distance [m]')
title('Fixed Frame')
legend('Y', 'Z')

subplot(3,1,2);
plot(t, wakecenterY_helixFrame_store)
hold on;
plot(t, wakecenterZ_helixFrame_store)
hold off;
xlabel('Time [s]')
ylabel('Distance [m]')
title('Helix Frame')
legend('Y', 'Z')

subplot(3,1,3);
plot(t, wakecenterY_fixFrame_store + mean(wakecenterZ_fixFrame_store))
hold on;
plot(t, wakecenterZ_fixFrame_store)
hold off;
xlabel('Time [s]')
ylabel('Distance [m]')
title('Again Fixed Frame')
legend('Y', 'Z')

%% Visualization
% ringVisualization(dataLiDAR)
