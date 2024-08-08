%% Helix of IEA15MW in script
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
fileName = 'FF_Uni_inflowAngle.mat';   % Fixed Frame
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 1000;     % in timestep, actual time is simTime*timestep(Q-blade define)
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
Turb_class = 'C';    % A, B, C
Turb_type = 'NTM';   % NTM, ETM, etc   
seed = 43;
vertInf = 0;         % Vertical inflow angle in degrees
horInf = 0;          % Horizontal inflow angle in degrees
% calllib('QBladeDLL', 'addTurbulentWind', ...
%     U_inflow,Hub_NREL5MW,Hub_NREL5MW,dimension,grid_point, ...
%     Turb_time,Turb_dt,Turb_class,Turb_type,seed,vertInf,horInf,1)

%% Defining Torque Control Setting
K = 2.24;
N = 97;          % Gearbox ratio

%% Defining Helix Control Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude = 3;                % Helix amplitude                
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;

% Genereate the tilt and yaw signal 
% t = linspace(timeStep, simLen, simTime);
t = linspace(1, simLen, simTime);
cutPoint = simLen/4;
cutPoint2 = simLen*3/4;
amplitudeTilt = Helix_amplitude*(t<=cutPoint) + 2*Helix_amplitude*(t>cutPoint) - Helix_amplitude*(t>=cutPoint2);
% plot(amplitudeTilt)
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
% legend('\beta_{tilt}','\beta_{yaw}')

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;   % 5(ring) to speed up sampling, only 4 valid points

%% Preparation for Simulation 
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
FF_theta = zeros(simTime, 2);
HF_theta = zeros(simTime, 2);
PitchAngles = zeros(simTime, 3);
FF_helixCenter = zeros(simTime, 2);
HF_helixCenter = zeros(simTime, 2);
templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
LiDAR_data(simTime, 1) = templateStruct;

% Sliding window
result = zeros(simTime, 2);
result_e = zeros(simTime, 2);
ws_filter = 100;
ws_centering = ceil(1/(Freq * timeStep));

% Deadtime Delay
timeDelay = LiDAR_x / U_inflow;

%% Real-time LPF
Fs = 1/timeStep;
Fc = 0.05;
Wn = Fc / (Fs / 2);

% Finite Impulse Response LPF (small phase lag in real-time)
n = 50; % Filter order
b_fir = fir1(n, Wn, 'low');
bufferSize = 50;   % 50
buffer = zeros(bufferSize, 2);
filterState1 = zeros(n, 1);
filterState2 = zeros(n, 1);
filterState3 = zeros(n, 1);
filterState4 = zeros(n, 1);

%% Simulation
% start simulation
tic
f = waitbar(0,'Initializing Simulation');
for i = 1:1:simTime
    calllib('QBladeDLL','advanceTurbineSimulation')
  
    % Get current value
    omega = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0.5, 0);
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
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);
    % 2. Inverse MBC 
    invMBC = [1 cosd(Azimuth1) sind(Azimuth1);
              1 cosd(Azimuth2) sind(Azimuth2);
              1 cosd(Azimuth3) sind(Azimuth3)];
    % 3. Blade pitch signal
    thetaBlade_Helix = invMBC * [0; 
                                 theta_tilt; 
                                 theta_yaw];   

    % 4. Corresponding helix frame value
    R_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
               -sin(omega_e*t(i)) cos(omega_e*t(i))];
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    thetaTiltYaw_helix = R_helix * [theta_tilt; 
                                    theta_yaw]; 

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        thetaBlade_Helix(1) thetaBlade_Helix(2) thetaBlade_Helix(3)],0)

    % LiDAR data sampling (Ring)
    windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample);
    wakeCenter = HelixCenter(windspeed, U_inflow, D_NREL5MW);
    FF_helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)

    % Get the helix center from the helix frame
    % LPF the single element
%     buffer = [buffer(2:end, :); FF_helixCenter(i, :)];
%     if i >= bufferSize
%         % Forward filtering
%         [forwardFiltered1, filterState1] = filter(b_fir, 1, buffer(:, 1), filterState1);
%         [forwardFiltered2, filterState2] = filter(b_fir, 1, buffer(:, 2), filterState2);
% 
%         % Reverse filtering
%         reverseFiltered1 = filter(b_fir, 1, flip(forwardFiltered1));
%         reverseFiltered2 = filter(b_fir, 1, flip(forwardFiltered2));
% 
%         % Reverse back to original order
%         result(i, 1) = reverseFiltered1(end);
%         result(i, 2) = reverseFiltered2(end);
%     else
%         % Initial phase where buffer is not yet full
%         result(i, 1) = FF_helixCenter(i, 1);
%         result(i, 2) = FF_helixCenter(i, 2);
%     end
    [result(i, 1), filterState1] = filter(b_fir, 1, FF_helixCenter(i, 1), filterState1);
    [result(i, 2), filterState2] = filter(b_fir, 1, FF_helixCenter(i, 2), filterState2);

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
%     center_e2 = invR_helix * [result(i, 1)-meanZ; result(i, 2)-meanY];
    % apply LPF to HF HelixCenter 
    [result_e(i, 1), filterState3] = filter(b_fir, 1, center_e(1), filterState3);
    [result_e(i, 2), filterState4] = filter(b_fir, 1, center_e(2), filterState4);

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
    LiDAR_data(i) = windspeed;

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
% plot(PitchAngles(:,1))
% hold on
% plot(PitchAngles(:,2))
% plot(PitchAngles(:,3))
% xticks(0:100:length(PitchAngles));
% xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
% xlabel("Time (s)");
% ylabel("Angle (deg)");
% legend('Blade 1','Blade 2','Blade 3')

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
plot(result(:, 1));
plot(result(:, 2));
hold off;
title('Center FF')
legend('z', 'y', 'z2', 'y2')
subplot(2, 2, 4)
plot(HF_helixCenter(:, 1));
hold on;
plot(HF_helixCenter(:, 2));
plot(result_e(:, 1));
plot(result_e(:, 2));
hold off;
title('Center HF')
legend('z_e', 'y_e', 'z_e2', 'y_e2')

% ringVisualization(LiDAR_data, D_NREL5MW)

%% Unload Library 
% unloadlibrary 'QBladeDLL'