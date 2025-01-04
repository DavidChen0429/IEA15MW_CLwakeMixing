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
simTime = 1680*2;     % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds
saveOption = 'Y';
HelixShape = 'basic'; % oval or flower

turbineName = '.\Data\NREL5MW\';
caseName = 'Sth\';
fileName = 'basic2.mat';
% QprName = 'oval.qpr';

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
Pit1 = 'Pitch Angle Blade 1 [deg]';
Pit2 = 'Pitch Angle Blade 2 [deg]';
Pit3 = 'Pitch Angle Blade 3 [deg]';
PowerVar = 'Aerodynamic Power [kW]';
CpVar = 'Power Coefficient [-]';
Moop1Var = 'Aero. OOP RootBend. Mom. Blade 1 [Nm]';
Mip1Var = 'Aero. IP RootBend. Mom. Blade 1 [Nm]';
Moop2Var = 'Aero. OOP RootBend. Mom. Blade 2 [Nm]';
Mip2Var = 'Aero. IP RootBend. Mom. Blade 2 [Nm]';
Moop3Var = 'Aero. OOP RootBend. Mom. Blade 3 [Nm]';
Mip3Var = 'Aero. IP RootBend. Mom. Blade 3 [Nm]';

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
AzimuthOffset = 96; % 6 for pi/2 shift ;96 for pi shift (right relationship)

% Genereate the tilt and yaw signal 
t = linspace(1, simLen, simTime);
t1 = linspace(1, simLen/4, simTime/4);
t2 = linspace(simLen/4, simLen/2, simTime/4);
t3 = linspace(simLen/2, simLen*3/4, simTime/4);
t4 = linspace(simLen*3/4, simLen, simTime/4);
if strcmp(HelixShape, 'oval')
    sigTilt = 3 * sin(2*pi*Freq*t);          
    sigYaw = 1 * sin(2*pi*Freq*t - pi/2);  % CCW
elseif strcmp(HelixShape, 'flower')
    sigTilt = [3*sin(2*pi*Freq*t1) 1*sin(2*pi*Freq*t2) 3*sin(2*pi*Freq*t3) 1*sin(2*pi*Freq*t4)];          
    sigYaw = [1*sin(2*pi*Freq*t1-pi/2) 3*sin(2*pi*Freq*t2-pi/2) 1*sin(2*pi*Freq*t3-pi/2) 3*sin(2*pi*Freq*t4-pi/2)];  % CCW
elseif strcmp(HelixShape, 'basic')
    sigTilt = 3 * sin(2*pi*Freq*t);          
    sigYaw = 3 * sin(2*pi*Freq*t - pi/2);  % CCW
end

Trigger = 0;

%% Defining LiDAR sampling 
% When you change this, don't forget to change the name of data.mat
LiDAR_x = 1*D_NREL5MW;   % Definition of x is pointing downwind
LiDAR_y = 0;
LiDAR_z = Hub_NREL5MW;   % Wind height
LiDAR_num_sample = 80;   % 5(ring) to speed up sampling, only 4 valid points

%% Preparation for Simulation 
% pre-define array to speed up code
TSR_store = zeros(simTime, 1);
Power_store = zeros(simTime, 1);
% blade 1
Moop1_store = zeros(simTime, 1);
Mip1_store = zeros(simTime, 1);
Mflap1_store = zeros(simTime, 1);
Medge1_store = zeros(simTime, 1);
% blade 2
Moop2_store = zeros(simTime, 1);
Mip2_store = zeros(simTime, 1);
Mflap2_store = zeros(simTime, 1);
Medge2_store = zeros(simTime, 1);
% blade 3
Moop3_store = zeros(simTime, 1);
Mip3_store = zeros(simTime, 1);
Mflap3_store = zeros(simTime, 1);
Medge3_store = zeros(simTime, 1);

Cp_store = zeros(simTime, 1);
FF_beta = zeros(simTime, 2);
HF_beta = zeros(simTime, 2);
FF_helixCenter_filtered = zeros(simTime, 2);
HF_helixCenter_filtered = zeros(simTime, 2);
PitchAngles = zeros(simTime, 3);
FF_helixCenter = zeros(simTime, 2);
HF_helixCenter = zeros(simTime, 2);
templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
LiDAR_data(simTime, 1) = templateStruct;

% Sliding window
ws_filter = 100;
ws_centering = ceil(1/(Freq * timeStep));

% Deadtime Delay
timeDelay = LiDAR_x / U_inflow;

%% Real-time LPF
[b_fir, n] = FIR_LPF(1/timeStep, 0.05);
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
    omega = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0, 0);
    genTorqueQB = calllib('QBladeDLL','getCustomData_at_num',valuestr2, 0, 0);
    TSR = calllib('QBladeDLL','getCustomData_at_num',valuestr3, 0, 0);
    Azimuth1 = calllib('QBladeDLL','getCustomData_at_num', Azimu1, 0, 0);
    Azimuth2 = calllib('QBladeDLL','getCustomData_at_num', Azimu2, 0, 0);
    Azimuth3 = calllib('QBladeDLL','getCustomData_at_num', Azimu3, 0, 0);
    Pitch1 = calllib('QBladeDLL','getCustomData_at_num', Pit1, 0, 0);
    Pitch2 = calllib('QBladeDLL','getCustomData_at_num', Pit2, 0, 0);
    Pitch3 = calllib('QBladeDLL','getCustomData_at_num', Pit3, 0, 0);
    Power = calllib('QBladeDLL','getCustomData_at_num', PowerVar, 0, 0);
    Cp = calllib('QBladeDLL','getCustomData_at_num', CpVar, 0, 0);
    Moop1 = calllib('QBladeDLL','getCustomData_at_num', Moop1Var, 0, 0);
    Mip1 = calllib('QBladeDLL','getCustomData_at_num', Mip1Var, 0, 0);
    Moop2 = calllib('QBladeDLL','getCustomData_at_num', Moop2Var, 0, 0);
    Mip2 = calllib('QBladeDLL','getCustomData_at_num', Mip2Var, 0, 0);
    Moop3 = calllib('QBladeDLL','getCustomData_at_num', Moop3Var, 0, 0);
    Mip3 = calllib('QBladeDLL','getCustomData_at_num', Mip3Var, 0, 0);
    
    % Define transform matrix 
    invMBC = [1 cosd(Azimuth1+AzimuthOffset) sind(Azimuth1+AzimuthOffset);
              1 cosd(Azimuth2+AzimuthOffset) sind(Azimuth2+AzimuthOffset);
              1 cosd(Azimuth3+AzimuthOffset) sind(Azimuth3+AzimuthOffset)];
    invR_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
                  -sin(omega_e*t(i)) cos(omega_e*t(i))];

    % ==================== LiDAR data sampling (Circle) 
    windspeed = Circle_LiDAR_Parallel(LiDAR_x, LiDAR_y, LiDAR_z, D_NREL5MW, LiDAR_num_sample); 
    wakeCenter = HelixCenter(windspeed, U_inflow, D_NREL5MW);
    FF_helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)
    % Get the helix center from the helix frame
    % LPF the single element
    [FF_helixCenter_filtered(i, 1), filterState1] = filter(b_fir, 1, FF_helixCenter(i, 1), filterState1);
    [FF_helixCenter_filtered(i, 2), filterState2] = filter(b_fir, 1, FF_helixCenter(i, 2), filterState2);

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
%     omega_store(i,:) = omega;
%     genTorqueQB_store(i,:) = genTorqueQB;
%     genTorque_store(i,:) = genTorque;
    TSR_store(i) = TSR;
    Power_store(i) = Power;
    Cp_store(i) = Cp;
    % blade 1
    Moop1_store(i) = Moop1;
    Mip1_store(i) = Mip1;
    Mflap1_store(i) = Moop1*cosd(Pitch1) + Mip1*sind(Pitch1);
    Medge1_store(i) = -Moop1*sind(Pitch1) + Mip1*cosd(Pitch1);
    % blade 2
    Moop2_store(i) = Moop2;
    Mip2_store(i) = Mip2;
    Mflap2_store(i) = Moop2*cosd(Pitch2) + Mip2*sind(Pitch2);
    Medge2_store(i) = -Moop2*sind(Pitch2) + Mip2*cosd(Pitch2);
    % blade 3
    Moop3_store(i) = Moop3;
    Mip3_store(i) = Mip3;
    Mflap3_store(i) = Moop3*cosd(Pitch3) + Mip3*sind(Pitch3);
    Medge3_store(i) = -Moop3*sind(Pitch3) + Mip3*cosd(Pitch3);
    
    FF_beta(i,:) = [theta_tilt theta_yaw];
%     AzimuthAngles(i,:) = [Azimuth1 Azimuth2 Azimuth3];
    PitchAngles(i,:) = [Pitch1 Pitch2 Pitch3];
    FF_helixCenter(i, :) = [wakeCenter(1) wakeCenter(2)]; % Z(tilt), Y(yaw)
    LiDAR_data(i) = windspeed;

    waitbar(i/simTime, f, sprintf('Simulation Running: %.1f%%', (i/simTime)*100));

end
close(f)
toc 
if strcmp(saveOption, 'Y')
%     calllib('QBladeDLL','storeProject', [turbineName caseName QprName]) 
%     save([turbineName caseName fileName], 'LiDAR_data', ...
%                                           'FF_helixCenter', ...
%                                           'FF_helixCenter_filtered', ...
%                                           'FF_beta', ...
%                                           'Power_store', ...
%                                           'Cp_store', ...
%                                           'Moop1_store', ...
%                                           'Mip1_store', ...
%                                           'Mflap1_store', ...
%                                           'Medge1_store', ...
%                                           'Moop2_store', ...
%                                           'Mip2_store', ...
%                                           'Mflap2_store', ...
%                                           'Medge2_store', ...
%                                           'Moop3_store', ...
%                                           'Mip3_store', ...
%                                           'Mflap3_store', ...
%                                           'Medge3_store', ...
%                                           'PitchAngles');
    save([turbineName caseName fileName], 'FF_helixCenter', ...
                                          'FF_helixCenter_filtered', ...
                                          'FF_beta', ...
                                          'PitchAngles');
end
calllib('QBladeDLL','closeInstance')

%% Visualization
trigger_time = Trigger * timeStep;

filter = 1000;
plot(FF_helixCenter_filtered(filter:end, 2), FF_helixCenter_filtered(filter:end, 1))
xlim([-24 14])
ylim([70 110])
% ringVisualization2(LiDAR_data, D_NREL5MW)