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
simFile = [SourcePath 'NREL5MW_1turbine.sim'];
addpath('.\Functions');

loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

%% Data file 
simTime = 3000;     % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds
saveOption = 'N';

turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\QBladeDeug\';
fileName = ['checkPhase_FF','.mat'];

%% Load project and Initialize simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')

% Variables we care
valuestr = 'Rotational Speed [rpm]';
valuestr2 = 'Gen. HSS Torque [Nm]';
valuestr3 = 'Tip Speed Ratio [-]';
Azimu1 = 'Azimuthal Position Blade 1 [deg]';
Azimu2 = 'Azimuthal Position Blade 2 [deg]';
Azimu3 = 'Azimuthal Position Blade 3 [deg]';

U_inflow = 10;        % Inflow wind speed, same with the Q-blade setting
D_NREL5MW = 126;     % Rotor diameter
Hub_NREL5MW = 90;   % Hub height

%% Defining Torque Control Setting
K = 2.24;
N = 97;          % Gearbox ratio

%% Defining Helix Control Setting
Str = 0.3;                          % Strouhal number
Helix_amplitude1 = 3;                % tilt 
Helix_amplitude2 = 0;                % yaw   
Freq = Str*U_inflow/D_NREL5MW;      % From Str, in Hz
omega_e = Freq*2*pi;
AzimuthOffset = 96;

% Genereate the tilt and yaw signal 
calibrate = 200;
t = linspace(calibrate, simLen+calibrate, simTime);
sigTilt = Helix_amplitude1 * sin(2*pi*Freq*t);          
sigYaw = Helix_amplitude2 * sin(2*pi*Freq*t - pi/2);  % CCW

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;   % 5(ring) to speed up sampling, only 4 valid points

%% Preparation for Simulation 
% T/4
FF_helixCenter1 = zeros(simTime, 2);
FF_helixCenter_filtered1 = zeros(simTime, 2);
HF_helixCenter_filtered1 = zeros(simTime, 2);
% 1D
FF_helixCenter2 = zeros(simTime, 2);
FF_helixCenter_filtered2 = zeros(simTime, 2);
HF_helixCenter_filtered2 = zeros(simTime, 2);
% T/2
FF_helixCenter3 = zeros(simTime, 2);
FF_helixCenter_filtered3 = zeros(simTime, 2);
HF_helixCenter_filtered3 = zeros(simTime, 2);

templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
LiDAR_data1(simTime, 1) = templateStruct;
LiDAR_data2(simTime, 1) = templateStruct;
LiDAR_data3(simTime, 1) = templateStruct;

ws_filter = 150;
ws_centering = ceil(1/(Freq * timeStep));

%% Real-time LPF
[b_fir, n] = FIR_LPF(1/timeStep, 0.05);
filterState1 = zeros(n, 1);
filterState2 = zeros(n, 1);
filterState3 = zeros(n, 1);
filterState4 = zeros(n, 1);
filterState5 = zeros(n, 1);
filterState6 = zeros(n, 1);

%% Simulation
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
    
    % Define transform matrix 
    invMBC = [1 cosd(Azimuth1+AzimuthOffset) sind(Azimuth1+AzimuthOffset);
              1 cosd(Azimuth2+AzimuthOffset) sind(Azimuth2+AzimuthOffset);
              1 cosd(Azimuth3+AzimuthOffset) sind(Azimuth3+AzimuthOffset)];
    invR_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
                  -sin(omega_e*t(i)) cos(omega_e*t(i))];

    % ==================== LiDAR data sampling (Circle) 
    windspeed1 = Circle_LiDAR_Parallel(105, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample); 
    wakeCenter1 = HelixCenter(windspeed1, U_inflow, D_NREL5MW);
    FF_helixCenter1(i, :) = [wakeCenter1(1) wakeCenter1(2)]; % Z(tilt), Y(yaw)
    [FF_helixCenter_filtered1(i, 1), filterState1] = filter(b_fir, 1, FF_helixCenter1(i, 1), filterState1);
    [FF_helixCenter_filtered1(i, 2), filterState2] = filter(b_fir, 1, FF_helixCenter1(i, 2), filterState2);

    windspeed2 = Circle_LiDAR_Parallel(126, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample); 
    wakeCenter2 = HelixCenter(windspeed2, U_inflow, D_NREL5MW);
    FF_helixCenter2(i, :) = [wakeCenter2(1) wakeCenter2(2)]; % Z(tilt), Y(yaw)
    [FF_helixCenter_filtered2(i, 1), filterState3] = filter(b_fir, 1, FF_helixCenter2(i, 1), filterState3);
    [FF_helixCenter_filtered2(i, 2), filterState4] = filter(b_fir, 1, FF_helixCenter2(i, 2), filterState4);

    windspeed3 = Circle_LiDAR_Parallel(210, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample); 
    wakeCenter3 = HelixCenter(windspeed3, U_inflow, D_NREL5MW);
    FF_helixCenter3(i, :) = [wakeCenter3(1) wakeCenter3(2)]; % Z(tilt), Y(yaw)
    [FF_helixCenter_filtered3(i, 1), filterState5] = filter(b_fir, 1, FF_helixCenter3(i, 1), filterState5);
    [FF_helixCenter_filtered3(i, 2), filterState6] = filter(b_fir, 1, FF_helixCenter3(i, 2), filterState6);
    
    meanZ = Hub_NREL5MW;
    meanY = 0;
    if i > ws_centering
        meanZ = mean(FF_helixCenter1(i-ws_centering:i, 1));
        meanY = mean(FF_helixCenter1(i-ws_centering:i, 2));
    end
    HF_helixCenter_filtered1(i, :) = invR_helix * [FF_helixCenter_filtered1(i, 1) - meanZ;
                                                   FF_helixCenter_filtered1(i, 2) - meanY];
    HF_helixCenter_filtered2(i, :) = invR_helix * [FF_helixCenter_filtered2(i, 1) - meanZ;
                                                   FF_helixCenter_filtered2(i, 2) - meanY];
    HF_helixCenter_filtered3(i, :) = invR_helix * [FF_helixCenter_filtered3(i, 1) - meanZ;
                                                   FF_helixCenter_filtered3(i, 2) - meanY];

    % ==================== Control
    % I. Torque control to maintain optimal TSR of 9 
    omega_g = omega*N;                      % rotor to generator
    genTorque = K.*(omega_g*(2*pi/60))^2;

    % II. Wake mixing
    % 1. Get tilt and yaw signals
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);
    % 2. Inverse MBC 
    % 3. Blade pitch signal
    thetaBlade_Helix = invMBC * [0; 
                                 theta_tilt; 
                                 theta_yaw];   

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        thetaBlade_Helix(1) thetaBlade_Helix(2) thetaBlade_Helix(3)],0)
    
    % ==================== Store values 
    LiDAR_data1(i) = windspeed1;
    LiDAR_data2(i) = windspeed2;
    LiDAR_data3(i) = windspeed3;
    waitbar(i/simTime, f, sprintf('Simulation Running: %.1f%%', (i/simTime)*100));

end
close(f)
toc 
if strcmp(saveOption, 'Y')
    save([turbineName caseName fileName], 'FF_helixCenter_filtered1', ...
                                          'FF_helixCenter_filtered2', ...
                                          'FF_helixCenter_filtered3');
end
calllib('QBladeDLL','closeInstance')

%% Visualization
filter = 1000;
figure('Position', [10, 10, 1200, 310]);
subplot(1, 3, 1)
plot(FF_helixCenter_filtered1(filter:end, 2), FF_helixCenter_filtered1(filter:end, 1));
xlim([-40 30])
ylim([57 127])
subplot(1, 3, 2)
plot(FF_helixCenter_filtered2(filter:end, 2), FF_helixCenter_filtered2(filter:end, 1));
xlim([-50 40])
ylim([47 137])
subplot(1, 3, 3)
plot(FF_helixCenter_filtered3(filter:end, 2), FF_helixCenter_filtered3(filter:end, 1));
xlim([-50 40])
ylim([47 137])

% 
% figure('Name', 'Overall Result', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% subplot(2, 1, 1)
% plot(HF_helixCenter_filtered1(filter:end, 1))
% hold on
% plot(HF_helixCenter_filtered2(filter:end, 1))
% plot(HF_helixCenter_filtered3(filter:end, 1))
% hold off
% title('z')
% legend('50','100','150')
% 
% subplot(2, 1, 2)
% plot(HF_helixCenter_filtered1(filter:end, 2))
% hold on
% plot(HF_helixCenter_filtered2(filter:end, 2))
% plot(HF_helixCenter_filtered3(filter:end, 2))
% hold off
% title('y')
% legend('50','100','150')

%% 
% figure('Position', [10, 450, 1200, 310]);
% interval = 10;
% for counter = 1:interval:simTime  
%     subplot(1, 3, 1)
%     snapshot = LiDAR_data1(counter);
%     y = snapshot.y;
%     z = snapshot.z;
%     scatter(y, z, 10, snapshot.u_los, 'filled');
%     xlabel('Y [m]')
%     ylabel('Z [m]')
%     title('LiDAR Wind Speed (sec)', round(counter/interval + 1))
%     colorbar;
%     clim([4 10])
% 
%     subplot(1, 3, 2)
%     snapshot2 = LiDAR_data2(counter);
%     y = snapshot2.y;
%     z = snapshot2.z;
%     scatter(y, z, 10, snapshot2.u_los, 'filled');
%     xlabel('Y [m]')
%     ylabel('Z [m]')
%     title('LiDAR Wind Speed (sec)', round(counter/interval + 1))
%     colorbar;
%     clim([4 10])
%     pause(0.1);
% 
%     subplot(1, 3, 3)
%     snapshot3 = LiDAR_data3(counter);
%     y = snapshot3.y;
%     z = snapshot3.z;
%     scatter(y, z, 10, snapshot3.u_los, 'filled');
%     xlabel('Y [m]')
%     ylabel('Z [m]')
%     title('LiDAR Wind Speed (sec)', round(counter/interval + 1))
%     colorbar;
%     clim([4 10])
%     pause(0.1);
% 
% end 