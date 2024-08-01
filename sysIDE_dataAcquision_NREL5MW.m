%% An attempt for LiDAR-enhanced-CL-Ctrl (Fixed Frame)
clear
close all 
%clc

%% Define paths
UserPath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\IEA15MW_CLwakeMixing\'; 
QBladePath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\QBladeEE_2.0.6.4\'; 
SourcePath = [UserPath 'Source\'];
DllPath = [QBladePath 'QBladeEE_2.0.6.dll'];
simFile = [SourcePath 'NREL5MW_Torque_Helix.sim'];
addpath('.\Functions');

loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

%% Data file 
fileName = 'biasZ1.mat';
dataPath = '.\Data\FF_bias\';
% caseName = 'FF_bias';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 10000;    % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds

% Variables we care
valuestr = 'Rotational Speed [rpm]';
valuestr2 = 'Gen. HSS Torque [Nm]';
valuestr3 = 'Tip Speed Ratio [-]';
Azimu1 = 'Azimuthal Position Blade 1 [deg]';
Azimu2 = 'Azimuthal Position Blade 2 [deg]';
Azimu3 = 'Azimuthal Position Blade 3 [deg]';
Pit1 = 'Pitch Angle Blade 1 [deg]';
Pit2 = 'Pitch Angle Blade 2 [deg]';
Pit3 = 'Pitch Angle Blade 3 [deg]';

%% Set Turbulent Wind
U_inflow = 10;        % Inflow wind speed, same with the Q-blade setting
D_NREL5MW = 126;     % Rotor diameter
Hub_NREL5MW = 90;   % Hub height
Wind_Height = Hub_NREL5MW;
dimension = D_NREL5MW;     % span dim*dim meters
grid_point = 50;     % sqaure grid
Turb_time = 10;      % Simulation length of the windfield in seconds
Turb_dt = timeStep;  % Temporal resolution of the windfield
Turb_class = 'A';    % A, B, C
Turb_type = 'NTM';   % NTM, ETM, etc   
seed = 43;
vertInf = 0;         % Vertical inflow angle in degrees
horInf = 0;          % Horizontal inflow angle in degrees
% calllib('QBladeDLL', 'addTurbulentWind', ...
%     U_inflow,Hub_IEA15MW,Hub_IEA15MW,dimension,grid_point, ...
%     Turb_time,Turb_dt,Turb_class,Turb_type,seed,vertInf,horInf,1)

%% Defining Torque Control Setting
% This need to be changed when inflow windspeed is varied
K = 2.24;
N = 97;          % Gearbox ratio

%% Defining Helix Control Setting
Str = 0.25;                          % Strouhal number
Helix_amplitude = 3;                % Helix amplitude                
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;

% Genereate the tilt and yaw signal 
t = linspace(timeStep, simLen, simTime);
% cutPoint = simLen/4;
% cutPoint2 = simLen*3/4;
% amplitudeTilt = Helix_amplitude*(t<=cutPoint) + 2*Helix_amplitude*(t>cutPoint) - Helix_amplitude*(t>=cutPoint2);

cutPoint = simLen/6;
sigTilt = Helix_amplitude * sin(2*pi*Freq*t);          
% sigTilt(t >= cutPoint) = sigTilt(t >= cutPoint) + 2;
sigYaw = Helix_amplitude * sin(2*pi*Freq*t - pi/2);  % CCW
% sigYaw(t >= cutPoint) = sigYaw(t >= cutPoint) - 5;

% sigTilt = amplitudeTilt .* sin(2*pi*Freq*t);  
% sigYaw = amplitudeTilt .* sin(2*pi*Freq*t - pi/2);  % CCW
% figure()
% plot(t, sigTilt)
% hold on
% plot(t, sigYaw)
% hold off
% xlabel('Time [s]')
% ylabel('Magnitude')
% legend('M_{tilt}','M_{yaw}')

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 60;   % 5(ring) to speed up sampling, only 4 valid points
DeadtimeDelay = LiDAR_x / U_inflow;

%% Signals for system IDE
n = 100; 
seedA = [1 0 0 1]; 
seedB = [1 1 0 1]; 
tap = [1 0 0 1]; % polynomial x^4 + x + 1

% PRBN Signal A
stateA = seedA;
prbnA = zeros(1, n);
for i = 1:n
    prbnA(i) = stateA(end); % Output the last bit
    feedbackA = mod(sum(stateA .* tap), 2); % Calculate feedback using tap positions
    stateA = [feedbackA stateA(1:end-1)]; % Shift and insert feedback
end

% PRBN Signal B
stateB = seedB;
prbnB = zeros(1, n);
for i = 1:n
    prbnB(i) = stateB(end); % Output the last bit
    feedbackB = mod(sum(stateB .* tap), 2); % Calculate feedback using tap positions
    stateB = [feedbackB stateB(1:end-1)]; % Shift and insert feedback
end

% scale 
% prbnA = prbnA * 5;
% prbnB = prbnB * 5;
prbnA(1:2) = 0;
prbnB(1:2) = 0;

% Plot the PRBN sequences
% figure;
% subplot(2, 1, 1);
% stairs(prbnA, 'LineWidth', 1.5);
% title('Pseudo-Random Binary Noise (PRBN) Sequence A');
% xlabel('Sample Index');
% ylabel('Binary Value');
% grid on;
% 
% subplot(2, 1, 2);
% stairs(prbnB, 'LineWidth', 1.5);
% title('Pseudo-Random Binary Noise (PRBN) Sequence B');
% xlabel('Sample Index');
% ylabel('Binary Value');
% grid on;

%% Ctrl Variables
% Ctrl inputs
% biasZArray = prbnA;
% biasYArray = prbnB;
% biasZArray = [0 0 0 1 1 1 1 1 1 1 1 1];
% biasYArray = [0 0 0 0 0 0 0 0 0 0 0 0];
biasZArray = [zeros(1, 5) 5*ones(1, 25)];
biasYArray = [zeros(1, 5) zeros(1, 25)];
biasZArray = randi([-5, 5], 1, 30);
biasZ = biasZArray(1);
biasY = biasYArray(1);
k = 2;

% Ctrl reference
RefY = 0;
RefZ = 90;

%% Simulation Data Array
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
thetaTiltYaw_fixedFrame = zeros(simTime, 2);
% thetaTiltYaw_helixFrame = zeros(simTime, 2);
PitchAngles = zeros(simTime, 3);
helixCenter = zeros(simTime, 2);
biasZStore = [];
biasYStore = [];
j = 0 ; % create for counting delay
% templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
% LiDAR_data(simTime, 1) = templateStruct;

% Sliding window to get the center of wake center
windowSize = 1/(Freq * timeStep);
buffer1 = [];
buffer2 = [];
meanValueY = [];
meanValueZ = [];

%% Start Simulation
tic
f = waitbar(0,'Initializing Simulation');
for i = 1:1:simTime
    calllib('QBladeDLL','advanceTurbineSimulation')
    j = j + 1;
    
    % Get current value
    omega = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0, 0);
    genTorqueQB = calllib('QBladeDLL','getCustomData_at_num',valuestr2, 0, 0);
    TSR = calllib('QBladeDLL','getCustomData_at_num',valuestr3, 0, 0);
    Azimuth1 = calllib('QBladeDLL','getCustomData_at_num', Azimu1, 0, 0);
    Azimuth2 = calllib('QBladeDLL','getCustomData_at_num', Azimu2, 0, 0);
    Azimuth3 = calllib('QBladeDLL','getCustomData_at_num', Azimu3, 0, 0);
    Pitch1 = calllib('QBladeDLL','getCustomData_at_num', Pit1, 0, 0);
    Pitch2 = calllib('QBladeDLL','getCustomData_at_num', Pit2, 0, 0);
    Pitch3 = calllib('QBladeDLL','getCustomData_at_num', Pit3, 0, 0);
    
    % Control
    % Torque control to maintain optimal TSR of 9 
    omega_g = omega*N;                      % rotor to generator
    genTorque = K.*(omega_g*(2*pi/60))^2;

    % Helix control
    % 1. Get tilt and yaw signals
    theta_col = 0;
    theta_tilt = sigTilt(i) + biasZ;
    theta_yaw = sigYaw(i) + biasY;
    % 2. Inverse MBC 
    invMBC = [1 cosd(Azimuth1) sind(Azimuth1);
              1 cosd(Azimuth2) sind(Azimuth2);
              1 cosd(Azimuth3) sind(Azimuth3)];
    % 3. Blade pitch signal
    thetaBlade_Helix = invMBC * [theta_col; 
                                 theta_tilt; 
                                 theta_yaw];   

    % 4. Corresponding helix frame value
%     R_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
%                -sin(omega_e*t(i)) cos(omega_e*t(i))];
%     thetaTiltYaw_helix = R_helix * [theta_tilt; theta_yaw]; 

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        thetaBlade_Helix(1) thetaBlade_Helix(2) thetaBlade_Helix(3)],0)

    % LiDAR data sampling   
    windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample); 
    wakeCenter = HelixCenter(windspeed, U_inflow, D_NREL5MW);

    % Store values 
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i) = TSR;
    thetaTiltYaw_fixedFrame(i,:) = [theta_tilt, theta_yaw];
%     thetaTiltYaw_helixFrame(i,:) = [thetaTiltYaw_helix(1), thetaTiltYaw_helix(2)];
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
    PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
%     LiDAR_data(i) = windspeed;
    helixCenter(i, :) = [wakeCenter(1), wakeCenter(2)];
%     bias1Store = [bias1Store bias1Array(k)];
%     bias2Store = [bias2Store bias2Array(k)];
    
    % calculating the center's center using sliding window
    if j >= DeadtimeDelay/timeStep
        buffer1 = [buffer1 wakeCenter(1)];
        buffer2 = [buffer2 wakeCenter(2)];
        if length(buffer1) == windowSize
            % Save wake center's center info
            buffer1 = lowpassFilter(buffer1, 1/timeStep, Freq+0.01);
            buffer2 = lowpassFilter(buffer2, 1/timeStep, Freq+0.01);
            meanValueY = [meanValueY mean(buffer1)];
            meanValueZ = [meanValueZ mean(buffer2)];
            biasZStore = [biasZStore biasZ];
            biasYStore = [biasYStore biasY];
           
            % Error signal
            errY = RefY - mean(buffer1);
            errZ = RefZ - mean(buffer2);
            errVec = [errZ errY];
            
            biasZ = biasZArray(k);  % 5 LTI?
            biasY = biasYArray(k);
            k = k + 1;
            
            % reset the sliding window
            buffer1 = [];  
            buffer2 = [];
            % if bias change, then recount the time delay
            if biasZ ~= biasZStore(end)
                j = 0;         
            end             
        end
    end

%     if mod(i, 1/timeStep) == 0
% %         fprintf('%d seconds.\n', i*timeStep);
%         LiDAR_data(i) = windspeed;   % 1Hz sampling
%     end
    waitbar(i/simTime, f, sprintf('Simulation Running: %.1f%%', (i/simTime)*100));

end
close(f)
%calllib('QBladeDLL','storeProject','15MW_Helix_Uni-U8_Str3.qpr') 
calllib('QBladeDLL','closeInstance')
% save([dataPath caseName fileName], 'LiDAR_data', ...
%                                    'thetaTiltYaw_fixedFrame', ...
%                                    'thetaTiltYaw_helixFrame');
toc 

%% Visualization
% figure;
% plot(TSR_store)
% xticks(0:100:length(TSR_store));
% xticklabels(0:100*timeStep:length(TSR_store)*timeStep);
% legend('TSR')
% xlabel("Time (s)");
% ylabel("TSR");

% figure;
% plot(genTorqueQB_store)
% hold on
% plot(genTorque_store)
% xticks(0:100:length(genTorqueQB_store));
% xticklabels(0:100*timeStep:length(genTorqueQB_store)*timeStep);
% xlabel("Time (s)");
% ylabel("Torque (Nm)")
% legend('QB HSS Torque','K omega^2')

figure;
subplot(2, 2, 1)
plot(thetaTiltYaw_fixedFrame(:, 1));
hold on;
plot(thetaTiltYaw_fixedFrame(:, 2));
hold off;
xticks(0:100:length(thetaTiltYaw_fixedFrame));
xlabel('Time (s)')
ylabel('Magnitude')
title('Fixed Frame')
legend('\theta_{tilt}', '\theta_{yaw}')

subplot(2, 2, 2)
plot(PitchAngles(:,1))
hold on
plot(PitchAngles(:,2))
plot(PitchAngles(:,3))
hold off
xticks(0:100:length(PitchAngles));
xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
xlabel("Time (s)");
ylabel("Angle (deg)");
legend('Blade 1','Blade 2','Blade 3')

subplot(2, 2, 3)
plot(helixCenter(:,1))
hold on
plot(helixCenter(:,2))
plot(lowpassFilter(helixCenter(:,1), 1/timeStep, Freq+0.01))
plot(lowpassFilter(helixCenter(:,2), 1/timeStep, Freq+0.01))
plot(biasZStore)
plot(biasYStore)
hold off
xticks(0:100:length(helixCenter));
xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
xlabel("Time (s)");
ylabel("Position (m)");
legend('Y','Z','Y_f','Z_f','biasZ','biasY')

subplot(2, 2, 4)
% plot(helixCenter(:,1), helixCenter(:,2))
plot(lowpassFilter(helixCenter(:,1), 1/timeStep, Freq+0.01), lowpassFilter(helixCenter(:,2), 1/timeStep, Freq+0.01))
hold on
plot(meanValueY(1:end), meanValueZ(1:end),'*')
hold off
xticks(0:100:length(helixCenter));
xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
xlim([-30 30])
ylim([60 120])
xlabel("Y (m)");
ylabel("Z (m)");

figure;
plot(meanValueZ)
hold on
plot(meanValueY)
plot(biasZStore)
plot(biasYStore)
hold off
legend('Z_m', 'Y_m', 'biasZ', 'biasY')
% save([dataPath fileName], 'meanValueZ', ...
%                           'meanValueY', ...                 
%                           'biasZStore', ...
%                           'biasYStore');

% ringVisualization(LiDAR_data, D_NREL5MW)

%% Unload Library 
% unloadlibrary 'QBladeDLL'