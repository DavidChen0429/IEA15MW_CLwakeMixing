%% An attempt for LiDAR-enhanced-CL-Ctrl (Fixed Frame)
clear
close all 
%clc

%% Define paths
UserPath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\IEA15MW_CLwakeMixing\'; 
QBladePath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\QBladeEE_2.0.6.4\'; 
SourcePath = [UserPath 'Source\'];
DllPath = [QBladePath 'QBladeEE_2.0.6.dll'];
simFile = [SourcePath 'IEA15MW_torque_Helix.sim'];
addpath('.\Functions');

% loadlibrary(DllPath,'.\QBladeDLLFunctions.h','alias','QBladeDLL') 
loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

%% Data file 
fileName = '600s_Center_FF_try_control.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')
simTime = 18000;     % in timestep, actual time is simTime*timestep(Q-blade define)
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
U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;     % Rotor diameter
Hub_IEA15MW = 150;   % Hub height
Wind_Height = Hub_IEA15MW;
dimension = D_IEA15MW;     % span dim*dim meters
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
K_list = [0.33358007e6 
          0.344363428e6 
          0.344363428e6 
          0.344363428e6]; % k accord with u_in = [7 8 9 10]
K = K_list(2);   % gain for k-omega^2, T/RS^2
N = 1;          % gearbox ratio of IEA15MW direct drive

%% Defining Helix Control Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude = 2;                % Helix amplitude                
Freq = Str*U_inflow/D_IEA15MW;      % From Str, in Hz
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
LiDAR_x = 1*D_IEA15MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Wind_Height;   % Wind height
LiDAR_num_sample = 60;   % 5(ring) to speed up sampling, only 4 valid points
DeadtimeDelay = LiDAR_x / U_inflow;

%% Simulation
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
thetaTiltYaw_fixedFrame = zeros(simTime, 2);
thetaTiltYaw_helixFrame = zeros(simTime, 2);
PitchAngles = zeros(simTime, 3);
helixCenter = zeros(simTime, 2);
bias = zeros(simTime, 2);
% templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
% LiDAR_data(simTime, 1) = templateStruct;

% Sliding window to get the center of wake center
windowSize = 1/(Freq * timeStep);
buffer1 = [];
buffer2 = [];
meanValueY = [];
meanValueZ = [];

% Ctrl inputs
bias1 = 0;
bias2 = 0;
% bias1Array = [1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
bias1Array = zeros(1, 100);
bias1Array(mod(1:100, 2) == 1) = 5;
k = 1;
% Ctrl reference
RefY = 0;
RefZ = 150;

% start simulation
tic
f = waitbar(0,'Initializing Simulation');
for i = 1:1:simTime
    calllib('QBladeDLL','advanceTurbineSimulation')
    
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
    omega_g = omega*N;  % rotor to generator
    genTorque = K.*(omega_g)^2;

    % Helix control
    % 1. Get tilt and yaw signals
    theta_col = 0;
    theta_tilt = sigTilt(i) + bias1;
    theta_yaw = sigYaw(i) + bias2;
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
    thetaTiltYaw_helix = R_helix * [theta_tilt; theta_yaw]; 

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        thetaBlade_Helix(1) thetaBlade_Helix(2) thetaBlade_Helix(3)],0)

    % LiDAR data sampling   
    windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, LiDAR_num_sample);    
    wakeCenter = HelixCenter(windspeed, U_inflow);

    % Store values 
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i) = TSR;
    thetaTiltYaw_fixedFrame(i,:) = [theta_tilt, theta_yaw];
    thetaTiltYaw_helixFrame(i,:) = [thetaTiltYaw_helix(1), thetaTiltYaw_helix(2)];
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
    PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
%     LiDAR_data(i) = windspeed;
    helixCenter(i, :) = [wakeCenter(1), wakeCenter(2)];
    bias(i, :) = [bias1, bias2];
    
    % calculating the center's center using sliding window
    if i >= DeadtimeDelay/timeStep
        buffer1 = [buffer1 wakeCenter(1)];
        buffer2 = [buffer2 wakeCenter(2)];
%         if length(buffer) > windowSize
%             buffer = buffer(end-windowSize+1:end);
%         end
        if length(buffer1) == windowSize
            buffer1 = lowpassFilter(buffer1, 1/timeStep, Freq+0.01);
            buffer2 = lowpassFilter(buffer2, 1/timeStep, Freq+0.01);
            meanValueY = [meanValueY mean(buffer1)];
            meanValueZ = [meanValueZ mean(buffer2)];
            
            % error signal
            errY = RefY - mean(buffer1);
            errZ = RefZ - mean(buffer2);
%             bias1 = errY/2;   % simple P controller
%             bias2 = errZ/2;   % simple P controller
%             bias1 = 0.5;   % simple P controller
%             bias2 = -2;   % simple P controller
            errVec = [errZ errY];
%             Mtilt_unit = [2.72-(-7.525) 172.66-158.54]/5;   % (2.049 2.824)
%             Myaw_unit = [2.72-(-19.82) 162.80-158.54]/5;    % (4.508 0.852)
%             B = [Mtilt_unit(:), Myaw_unit(:)];
%             coefficients = B \ errVec(:);
%             bias1 = bias1 + coefficients(1);
%             bias2 = bias2 + coefficients(2);
            bias1 = bias1Array(k);
            bias2 = bias1Array(k);
            k = k + 1;
%             disp(errVec)
%             disp([bias1 bias2])
            
            buffer1 = [];   % reset the sliding window
            buffer2 = [];
%             buffer = buffer(overlap+1 : end);   % reset the sliding window
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

figure;
plot(helixCenter(:,1))
hold on
plot(helixCenter(:,2))
plot(lowpassFilter(helixCenter(:,1), 1/timeStep, Freq+0.01))
plot(lowpassFilter(helixCenter(:,2), 1/timeStep, Freq+0.01))
plot(bias(:, 1))
plot(bias(:, 2))
hold off
xticks(0:100:length(helixCenter));
xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
xlabel("Time (s)");
ylabel("Angle (deg)");
legend('Y','Z','Y_f','Z_f','bias1','bias2')

figure;
% plot(helixCenter(:,1), helixCenter(:,2))
plot(lowpassFilter(helixCenter(:,1), 1/timeStep, Freq+0.01), lowpassFilter(helixCenter(:,2), 1/timeStep, Freq+0.01))
hold on
% plot(lowpassFilter(helixCenter(:,1), 1/timeStep, Freq+0.01))
% plot(lowpassFilter(helixCenter(:,2), 1/timeStep, Freq+0.01))
plot(meanValueY, meanValueZ,'*')
hold off
xticks(0:100:length(helixCenter));
xticklabels(0:100*timeStep:length(PitchAngles)*timeStep);
xlim([-50 50])
ylim([100 200])
xlabel("Y (m)");
ylabel("Z (m)");
% legend('Y','Z','Y_f','Z_f')

% figure;
% plot(AzimuthAngles(:,1))
% hold on
% plot(AzimuthAngles(:,2))
% plot(AzimuthAngles(:,3))
% hold off
% xticks(0:100:length(AzimuthAngles));
% xticklabels(0:100*timeStep:length(AzimuthAngles)*timeStep);
% xlabel("Time (s)");
% ylabel("Angle (deg)");
% legend('Blade 1','Blade 2','Blade 3')

% figure;
% subplot(2, 1, 1)
% plot(thetaTiltYaw_fixedFrame(:, 1));
% hold on;
% plot(thetaTiltYaw_fixedFrame(:, 2));
% hold off;
% title('Fixed Frame')
% legend('\theta_{tilt}', '\theta_{yaw}')
% subplot(2, 1, 2);
% plot(thetaTiltYaw_helixFrame(:, 1));
% hold on;
% plot(thetaTiltYaw_helixFrame(: ,2));
% hold off;
% title('Helix Frame')
% legend('\theta^e_{tilt}', '\theta^e_{yaw}')

% ringVisualization(LiDAR_data)

%% Unload Library 
% unloadlibrary 'QBladeDLL'