%% An attempt for LiDAR-enhanced-CL-Ctrl (Fixed Frame)
clear
close all 
addpath('.\Functions');
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
fileName = 'FF_Uni_tilt,bias.mat';   % Fixed Frame
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 6000;    % in timestep, actual time is simTime*timestep(Q-blade define)
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
Str = 0.3;                          % Strouhal number
Helix_amplitude = 3;                % Helix amplitude                
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;

% Genereate the tilt and yaw signal 
t = linspace(1, simLen, simTime);
% t = linspace(timeStep, simLen, simTime);
sigTilt = Helix_amplitude * sin(2*pi*Freq*t);          
sigYaw = Helix_amplitude * sin(2*pi*Freq*t - pi/2);  % CCW

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;   % 5(ring) to speed up sampling, only 4 valid points
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
% biasZArray = [zeros(1, 2) 0*ones(1, 25)];
% biasYArray = [zeros(1, 2) -3*ones(1, 25)];
biasZArray = [zeros(1, 2) randi([-3, 3], 1, 30)];
biasYArray = [zeros(1, 2) randi([-3, 3], 1, 30)];
biasZ = biasZArray(1);
biasY = biasYArray(1);
k = 2;

% Ctrl reference
% RefY = 0;
% RefZ = 90;

%% Simulation Data Array
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
FF_theta = zeros(simTime, 2);
HF_theta = zeros(simTime, 2);
PitchAngles = zeros(simTime, 3);
FF_helixCenter = zeros(simTime, 2);
HF_helixCenter = zeros(simTime, 2);
biasZStore = [];
biasYStore = [];
j = 0 ; % create for counting delay
templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
LiDAR_data(simTime, 1) = templateStruct;

% Sliding window to get the center of wake center
buffer1 = [];
buffer2 = [];
meanValueY = [];
meanValueZ = [];

% Sliding window
result = zeros(simTime, 2);
result_e = zeros(simTime, 2);
ws_filter = 100;
ws_centering = ceil(1/(Freq * timeStep));

% Deadtime Delay
timeDelay = LiDAR_x / U_inflow;

% Low pass filter property
Fs = 1/timeStep;
Fc = 0.03;

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
    R_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
               -sin(omega_e*t(i)) cos(omega_e*t(i))];
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    thetaTiltYaw_helix = R_helix * [theta_tilt - biasZ; 
                                    theta_yaw - biasY]; 

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        thetaBlade_Helix(1) thetaBlade_Helix(2) thetaBlade_Helix(3)],0)

    % LiDAR data sampling   
    windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample); 
    wakeCenter = HelixCenter(windspeed, U_inflow, D_NREL5MW);
    FF_helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)
    
    % Get the mean
    meanZ = Hub_NREL5MW;
    meanY = 0;
    if i > ws_centering
        meanZ = mean(FF_helixCenter(i-ws_centering:i, 1));
        meanY = mean(FF_helixCenter(i-ws_centering:i, 2));
    end
    
    % Low pass filter
    % Centering
    centerZ = wakeCenter(1) - meanZ;  % 91.9411
    centerY = wakeCenter(2) - meanY;  % -3.1245
    center_e = invR_helix * [centerZ; centerY];

    % Store values 
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i) = TSR;
    FF_theta(i,:) = [theta_tilt theta_yaw];
    HF_theta(i,:) = [thetaTiltYaw_helix(1) thetaTiltYaw_helix(2)];
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
    PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
    FF_helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)
    HF_helixCenter(i, :) = [center_e(1) center_e(2)];   % Ze(tilt), Ye(yaw) 
%     result_e(i, :) = [center_e2(1) center_e2(2)];
    LiDAR_data(i) = windspeed;
    
    % calculating the center's center using sliding window
    if j >= DeadtimeDelay/timeStep
        buffer1 = [buffer1 wakeCenter(1)];
        buffer2 = [buffer2 wakeCenter(2)];
        if length(buffer1) == ws_centering
            % Save wake center's center info
            buffer1 = lowpassFilter(buffer1, 1/timeStep, Freq+0.01);
            buffer2 = lowpassFilter(buffer2, 1/timeStep, Freq+0.01);
            meanValueY = [meanValueY mean(buffer2)];
            meanValueZ = [meanValueZ mean(buffer1)];
            biasZStore = [biasZStore biasZ];
            biasYStore = [biasYStore biasY];
           
%             % Error signal
%             errY = RefY - mean(buffer1);
%             errZ = RefZ - mean(buffer2);
%             errVec = [errZ errY];
            
            biasZ = biasZArray(k);
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
% save([turbineName caseName fileName], 'LiDAR_data', ...
%                                       'FF_helixCenter', ...
%                                       'HF_helixCenter', ...
%                                       'FF_theta', ...
%                                       'HF_theta');
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
plot(FF_theta(:, 1));
hold on;
plot(FF_theta(:, 2));
hold off;
title('\beta FF')
legend('\theta_{tilt}', '\theta_{yaw}')
subplot(2, 2, 3);
plot(HF_theta(:, 1));
hold on;
plot(HF_theta(: ,2));
hold off;
title('\beta_e HF')
legend('\theta^e_{tilt}', '\theta^e_{yaw}')
subplot(2, 2, 2)
plot(FF_helixCenter(:, 1));
hold on;
plot(FF_helixCenter(:, 2));
hold off;
title('Center FF')
legend('z', 'y')
subplot(2, 2, 4)
% plot(HF_helixCenter(:, 1));
% hold on;
% plot(HF_helixCenter(:, 2));
plot(lowpassFilter(HF_helixCenter(:, 1),1/timeStep,Freq+0.01));
hold on;
plot(lowpassFilter(HF_helixCenter(:, 2),1/timeStep,Freq+0.01));
hold off;
title('Center HF')
legend('z_e', 'y_e')
subplot(2, 2, 4)

figure('Position', [10, 10, 500, 500]);
plot(lowpassFilter(FF_helixCenter(:,2), 1/timeStep, Freq+0.01), lowpassFilter(FF_helixCenter(:,1), 1/timeStep, Freq+0.01))
hold on
plot(meanValueY(1:end), meanValueZ(1:end),'*')
hold off
xlim([-30 30])
ylim([60 120])
xlabel("Y (m)");
ylabel("Z (m)");

% ringVisualization(LiDAR_data, D_NREL5MW)

%% Unload Library 
% unloadlibrary 'QBladeDLL'