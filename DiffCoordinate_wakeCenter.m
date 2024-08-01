%% This file is created to learn the different cooridnate for the helix center
clear
close all
addpath('.\Functions');

%% Basic information definition
fileName = '600s_Center_FF_baseline.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';
SimData = load([dataPath caseName fileName]);
% wakeCenterChange_Visualization(SimData)
dataLiDAR= SimData.LiDAR_data;
filteredIndex = 1;
dataLiDAR = dataLiDAR(filteredIndex:end);
data_length = size(dataLiDAR);              % length of snapshot
lengthPoint = length(dataLiDAR(1).x);

U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;     % Rotor diameter
simTime = data_length(1);     % in timestep, actual time is simTime*timestep(Q-blade define)
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

%% Perfect wake center
idealZ = 25 * ones(1, data_length(1));  % 25
idealY = -10 * ones(1, data_length(1)); % -10

% idealZ = [25*ones(1, data_length(1)/4) 50*ones(1, data_length(1)/2) 25*ones(1, data_length(1)/4)];
% idealY = [-10*ones(1, data_length(1)/4) 0*ones(1, data_length(1)/2) -10*ones(1, data_length(1)/4)];
% idealZ = linspace(0, 50, data_length(1));

%% From Fixed Frame --> Helix Frame, CCW
% Fixed Frame
wakeCenterY = arrayfun(@(x) x.centerY, dataLiDAR);
wakeCenterZ = arrayfun(@(x) x.centerZ, dataLiDAR);

% LPF
wakeCenterY_f = lowpassFilter(wakeCenterY, 10, 0.02);
wakeCenterZ_f = lowpassFilter(wakeCenterZ, 10, 0.02);
% figure()
% plot(wakeCenterY_f)
% hold on
% plot(wakeCenterZ_f)
% hold off

% Normalize to get rid of previous 1D gain
wakeCenterY = wakeCenterY_f - mean(wakeCenterY_f);
wakeCenterZ = wakeCenterZ_f - mean(wakeCenterZ_f);
t = linspace(1, simLen, data_length(1));

% Store wake center information in helix frame
wakecenterY_helixFrame_store = zeros(data_length(1), 1);
wakecenterZ_helixFrame_store = zeros(data_length(1), 1);
wakecenterY_fixFrame_store = zeros(data_length(1), 1);
wakecenterZ_fixFrame_store = zeros(data_length(1), 1);
ideal_wakecenterY_fixFrame_store = zeros(data_length(1), 1);
ideal_wakecenterZ_fixFrame_store = zeros(data_length(1), 1);

for i = 1:1:data_length(1)
    R_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
               sin(omega_e*t(i)) cos(omega_e*t(i))];
    invR_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
                  -sin(omega_e*t(i)) cos(omega_e*t(i))];
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

    ideal_center_fixedFrame = invR_helix * [idealZ(i); idealY(i)];
    ideal_wakecenterY_fixFrame_store(i) = ideal_center_fixedFrame(2);
    ideal_wakecenterZ_fixFrame_store(i) = ideal_center_fixedFrame(1);
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
plot(t, idealY, '--')
plot(t, idealZ, '--')
hold off;
xlabel('Time [s]')
ylabel('Distance [m]')
title('Helix Frame')
legend('Y^e', 'Z^e', 'Y^e_r', 'Z^e_r')

subplot(3,1,3);
plot(t, wakecenterY_fixFrame_store + mean(wakecenterZ_fixFrame_store))
hold on;
plot(t, wakecenterZ_fixFrame_store)
plot(t, ideal_wakecenterY_fixFrame_store, '--')
plot(t, ideal_wakecenterZ_fixFrame_store, '--')
hold off;
xlabel('Time [s]')
ylabel('Distance [m]')
title('Again Fixed Frame')
legend('Y', 'Z', 'Y_r', 'Z_r')

%% Visualization
figure('Position', [10, 10, 500, 500]);
plot(wakeCenterY, wakeCenterZ,'red');
hold on
plot(mean(wakeCenterY), mean(wakeCenterZ),'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(ideal_wakecenterY_fixFrame_store, ideal_wakecenterZ_fixFrame_store,'blue')
plot(mean(ideal_wakecenterY_fixFrame_store), mean(ideal_wakecenterZ_fixFrame_store),'bo', 'MarkerSize', 10, 'LineWidth', 2);
hold off
xlabel('Y [m]')
ylabel('Z [m]')
xlim([-50 50])
ylim([-50 50])
title('Wake center Trajectory')