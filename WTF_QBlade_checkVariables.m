%% Visualize Metrices for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Define paths
UserPath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\IEA15MW_CLwakeMixing\'; 
QBladePath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\QBladeEE_2.0.6.4\'; 
SourcePath = [UserPath 'Source\'];
DllPath = [QBladePath 'QBladeEE_2.0.6.dll'];
simFile = [SourcePath 'NREL5MW_2turbines_4D.sim'];
addpath('.\Functions');

loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

%% Basic Simulation Settings
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\QBladeDeug\';

% Simulation Settings
simTime = 6000;    
timeStep = 0.1;    
simLen = simTime * timeStep; 
windCase = 'Uniform';   % Basline, Uniform
controlStrategy = 'OL'; % OL CL
saveOption = 'Y';

% Turbine Settings
U_inflow = 10;
D_NREL5MW = 126;
Hub_NREL5MW = 90;

% LiDAR Settings
MeasurementPos = 4*D_NREL5MW-50;

% Control Settings
K = 2.24;        % NREL5MW
N = 97;          % Gearbox ratio

% File Settings
InputfileName = 'inputList.mat';
savefileName = ['2Turbines_',windCase,'_',controlStrategy,'.mat'];

%% Load Input Data
if strcmp(windCase, 'Baseline')
    InputSignal.Pitch = zeros(simTime, 3);
else
    InputBuf = load([turbineName caseName InputfileName]);
    InputSignal = InputBuf.InputData.(windCase).(controlStrategy);
end

%% Simulation
%this is setup using relative path and depends on the location of this file
calllib('QBladeDLL','createInstance',2,64)  % 64 for ring
calllib('QBladeDLL','setLibraryPath',DllPath)   % set lib path
calllib('QBladeDLL','loadSimDefinition',simFile)
calllib('QBladeDLL','initializeSimulation')

% Variables we care
valuestr = 'Rotational Speed [rpm]';
valuestr2 = 'Gen. HSS Torque [Nm]';
valuestr3 = 'Tip Speed Ratio [-]';

TorqueStoreTurb1 = zeros(simTime, 1);
TorqueStoreTurb2 = zeros(simTime, 1);
templateStruct = struct('x', [], 'y', [], 'z', [], 'u_x', [], 'u_y', [], 'u_z', [], 'u_norm', [], 'u_los', []);
WindData(simTime, 1) = templateStruct;
HelixTest(simTime, 1) = templateStruct;

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

    omega_turb2 = calllib('QBladeDLL','getCustomData_at_num',valuestr, 0, 1);
    genTorqueQB_turb2 = calllib('QBladeDLL','getCustomData_at_num',valuestr2, 0, 1);
    TSR_turb2 = calllib('QBladeDLL','getCustomData_at_num',valuestr3, 0, 1);

    % ==================== LiDAR data sampling (Circle) 
    windspeed = Circle_LiDAR_Parallel(MeasurementPos, 0, Hub_NREL5MW, D_NREL5MW, 80); 
    helixWind = Circle_LiDAR_Parallel(D_NREL5MW, 0, Hub_NREL5MW, D_NREL5MW, 80); 
    
    % ==================== Control
    % I. Torque control to maintain optimal TSR of 9 
    omega_g = omega*N;                      % rotor to generator
    genTorque = K.*(omega_g*(2*pi/60))^2;
    omega_g_turb2 = omega_turb2*N;          % rotor to generator
    genTorque_turb2 = K.*(omega_g_turb2*(2*pi/60))^2;

    % Send control signal to qblade
    calllib('QBladeDLL','setControlVars_at_num',[genTorque 0 ...
        InputSignal.Pitch(i, 1) InputSignal.Pitch(i, 2) InputSignal.Pitch(i, 3)],0)
    calllib('QBladeDLL','setControlVars_at_num',[genTorque_turb2 0 ...
        0 0 0],1)

    % ==================== Store values 
    TorqueStoreTurb1(i, :) = genTorqueQB;
    TorqueStoreTurb2(i, :) = genTorqueQB_turb2;
    WindData(i) = windspeed;
    HelixTest(i) = helixWind;

    waitbar(i/simTime, f, sprintf('Simulation Running: %.1f%%', (i/simTime)*100));

end
close(f)
if strcmp(saveOption, 'Y')
    save([turbineName caseName savefileName], "TorqueStoreTurb1", ...
                                              "TorqueStoreTurb2", ...
                                              "WindData");
end
calllib('QBladeDLL','closeInstance')
toc 

%% See Data
% videoCompare_func2(HelixTest, WindData, D_NREL5MW, '');
