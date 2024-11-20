%% Visualize Metrices for Experiments
clear
close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2Turbines\';

% Different case
windCase = 'Uniform'; % Uniform, Shear, Turb, Both, Shear2
if strcmp(windCase, 'Uniform')
    basefile = '2Turbines_Baseline.mat';
    OLfileName = '2Turbines_OL_Helix_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_ramp&stop_mag3.mat';
elseif strcmp(windCase, 'Shear')
    basefile = '2Turbines_Baseline_Shear0.2.mat';
    OLfileName = '2Turbines_OL_Helix_Shear0.2_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_Shear0.2_mag3.mat';
elseif strcmp(windCase, 'Turb')
    basefile = '2Turbines_Baseline_TI6.mat';
    OLfileName = '2Turbines_OL_Helix_TI6_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_TI6_mag3.mat';    
elseif strcmp(windCase, 'Both')
    basefile = '2Turbines_Baseline_TI6&Shear0.2.mat';
    OLfileName = '2Turbines_OL_Helix_TI6&Shear0.2_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_TI6&Shear0.2_mag3.mat';   
elseif strcmp(windCase, 'Shear2')
    basefile = '2Turbines_Baseline_Shear0.2Cheat.mat';
    OLfileName = '2Turbines_OL_Helix_Shear0.2_mag3.mat';
    CLfileName = '2Turbines_CL_Helix_Shear0.2_mag3.mat';    
end

Baseline = load([turbineName caseName basefile]);
OL = load([turbineName caseName OLfileName]);
CL = load([turbineName caseName CLfileName]);

%% Define Basic Helix Information
D_NREL5MW = 126;
U_inflow = 10;
timeStep = 0.1;
filter = 3500;  % Steady-state value
simLength = length(Baseline.Cp_store);

%% Value Calculation
% ============== Baseline
% Upstream WT1
BL_result.WT1.power = calculatePower(filter,Baseline.Cp_store,D_NREL5MW,U_inflow); % [MW]
BL_result.WT1.DEL = calculateDEL(filter, ...
    Baseline.Mflap1_store,Baseline.Medge1_store, ...
    Baseline.Mflap2_store,Baseline.Medge2_store, ...
    Baseline.Mflap3_store,Baseline.Medge3_store, ...
    timeStep); % [Nm]
BL_result.WT1.PBD = calculatePBD(filter,Baseline.PitchAngles, ...
    Baseline.Mflap1_store,Baseline.Medge1_store, ...
    Baseline.Mflap2_store,Baseline.Medge2_store, ...
    Baseline.Mflap3_store,Baseline.Medge3_store); % [kNm deg]
% Downstream WT2
BL_result.WT2.power = calculatePower(filter,Baseline.Cpturb2_store,D_NREL5MW,U_inflow); % [MW]
BL_result.WT2.DEL = calculateDEL(filter, ...
    Baseline.Mflap1turb2_store,Baseline.Medge1turb2_store, ...
    Baseline.Mflap2turb2_store,Baseline.Medge2turb2_store, ...
    Baseline.Mflap3turb2_store,Baseline.Medge3turb2_store, ...
    timeStep); % [Nm]
BL_result.WT2.PBD = calculatePBD(filter,Baseline.PitchAnglesturb2, ...
    Baseline.Mflap1turb2_store,Baseline.Medge1turb2_store, ...
    Baseline.Mflap2turb2_store,Baseline.Medge2turb2_store, ...
    Baseline.Mflap3turb2_store,Baseline.Medge3turb2_store); % [kNm deg]
BL_result.WT1.DEL_flapwise = BL_result.WT1.DEL(1)+BL_result.WT1.DEL(3)+BL_result.WT1.DEL(5);
BL_result.WT1.DEL_edgewise = BL_result.WT1.DEL(2)+BL_result.WT1.DEL(4)+BL_result.WT1.DEL(6);
BL_result.WT2.DEL_flapwise = BL_result.WT2.DEL(1)+BL_result.WT2.DEL(3)+BL_result.WT2.DEL(5);
BL_result.WT2.DEL_edgewise = BL_result.WT2.DEL(2)+BL_result.WT2.DEL(4)+BL_result.WT2.DEL(6);
BL_result.All.power = BL_result.WT1.power + BL_result.WT2.power;
BL_result.All.DEL_flapwise = BL_result.WT1.DEL_flapwise + BL_result.WT2.DEL_flapwise;
BL_result.All.DEL_edgewise = BL_result.WT1.DEL_edgewise + BL_result.WT2.DEL_edgewise;
BL_result.All.PBD = BL_result.WT1.PBD + BL_result.WT2.PBD;

% ============== Open-Loop
OL_result.WT1.power = calculatePower(filter,OL.Cp_store,D_NREL5MW,U_inflow); % [MW]
OL_result.WT1.DEL = calculateDEL(filter, ...
    OL.Mflap1_store,OL.Medge1_store, ...
    OL.Mflap2_store,OL.Medge2_store, ...
    OL.Mflap3_store,OL.Medge3_store, ...
    timeStep); % [Nm]
OL_result.WT1.PBD = calculatePBD(filter,OL.PitchAngles, ...
    OL.Mflap1_store,OL.Medge1_store, ...
    OL.Mflap2_store,OL.Medge2_store, ...
    OL.Mflap3_store,OL.Medge3_store); % [kNm deg]
% Downstream WT2
OL_result.WT2.power = calculatePower(filter,OL.Cpturb2_store,D_NREL5MW,U_inflow); % [MW]
OL_result.WT2.DEL = calculateDEL(filter, ...
    OL.Mflap1turb2_store,OL.Medge1turb2_store, ...
    OL.Mflap2turb2_store,OL.Medge2turb2_store, ...
    OL.Mflap3turb2_store,OL.Medge3turb2_store, ...
    timeStep); % [Nm]
OL_result.WT2.PBD = calculatePBD(filter,OL.PitchAnglesturb2, ...
    OL.Mflap1turb2_store,OL.Medge1turb2_store, ...
    OL.Mflap2turb2_store,OL.Medge2turb2_store, ...
    OL.Mflap3turb2_store,OL.Medge3turb2_store); % [kNm deg]
OL_result.WT1.DEL_flapwise = OL_result.WT1.DEL(1)+OL_result.WT1.DEL(3)+OL_result.WT1.DEL(5);
OL_result.WT1.DEL_edgewise = OL_result.WT1.DEL(2)+OL_result.WT1.DEL(4)+OL_result.WT1.DEL(6);
OL_result.WT2.DEL_flapwise = OL_result.WT2.DEL(1)+OL_result.WT2.DEL(3)+OL_result.WT2.DEL(5);
OL_result.WT2.DEL_edgewise = OL_result.WT2.DEL(2)+OL_result.WT2.DEL(4)+OL_result.WT2.DEL(6);
OL_result.All.power = OL_result.WT1.power + OL_result.WT2.power;
OL_result.All.DEL_flapwise = OL_result.WT1.DEL_flapwise + OL_result.WT2.DEL_flapwise;
OL_result.All.DEL_edgewise = OL_result.WT1.DEL_edgewise + OL_result.WT2.DEL_edgewise;
OL_result.All.PBD = OL_result.WT1.PBD + OL_result.WT2.PBD;

% ============== Closed-Loop
CL_result.WT1.power = calculatePower(filter,CL.Cp_store,D_NREL5MW,U_inflow); % [MW]
CL_result.WT1.DEL = calculateDEL(filter, ...
    CL.Mflap1_store,CL.Medge1_store, ...
    CL.Mflap2_store,CL.Medge2_store, ...
    CL.Mflap3_store,CL.Medge3_store, ...
    timeStep); % [Nm]
CL_result.WT1.PBD = calculatePBD(filter,CL.PitchAngles, ...
    CL.Mflap1_store,CL.Medge1_store, ...
    CL.Mflap2_store,CL.Medge2_store, ...
    CL.Mflap3_store,CL.Medge3_store); % [kNm deg]
% Downstream WT2
CL_result.WT2.power = calculatePower(filter,CL.Cpturb2_store,D_NREL5MW,U_inflow); % [MW]
CL_result.WT2.DEL = calculateDEL(filter, ...
    CL.Mflap1turb2_store,CL.Medge1turb2_store, ...
    CL.Mflap2turb2_store,CL.Medge2turb2_store, ...
    CL.Mflap3turb2_store,CL.Medge3turb2_store, ...
    timeStep); % [Nm]
CL_result.WT2.PBD = calculatePBD(filter,CL.PitchAnglesturb2, ...
    CL.Mflap1turb2_store,CL.Medge1turb2_store, ...
    CL.Mflap2turb2_store,CL.Medge2turb2_store, ...
    CL.Mflap3turb2_store,CL.Medge3turb2_store); % [kNm deg]
CL_result.WT1.DEL_flapwise = CL_result.WT1.DEL(1)+CL_result.WT1.DEL(3)+CL_result.WT1.DEL(5);
CL_result.WT1.DEL_edgewise = CL_result.WT1.DEL(2)+CL_result.WT1.DEL(4)+CL_result.WT1.DEL(6);
CL_result.WT2.DEL_flapwise = CL_result.WT2.DEL(1)+CL_result.WT2.DEL(3)+CL_result.WT2.DEL(5);
CL_result.WT2.DEL_edgewise = CL_result.WT2.DEL(2)+CL_result.WT2.DEL(4)+CL_result.WT2.DEL(6);
CL_result.All.power = CL_result.WT1.power + CL_result.WT2.power;
CL_result.All.DEL_flapwise = CL_result.WT1.DEL_flapwise + CL_result.WT2.DEL_flapwise;
CL_result.All.DEL_edgewise = CL_result.WT1.DEL_edgewise + CL_result.WT2.DEL_edgewise;
CL_result.All.PBD = CL_result.WT1.PBD + CL_result.WT2.PBD;

%% Visualize result
lw = 2;
Font = 20;

% === Power 
x = [1 2 3];
deltaPower = [(OL_result.WT1.power-BL_result.WT1.power)/(BL_result.WT1.power) (OL_result.WT2.power-BL_result.WT2.power)/(BL_result.WT2.power) (OL_result.All.power-BL_result.All.power)/(BL_result.All.power);
              (CL_result.WT1.power-BL_result.WT1.power)/(BL_result.WT1.power) (CL_result.WT2.power-BL_result.WT2.power)/(BL_result.WT2.power) (CL_result.All.power-BL_result.All.power)/(BL_result.All.power)];
figure('Name', 'Power', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
bar(x,deltaPower);
xticks(x); 
xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
ylabel('\Delta Power [%]')
legend('OL', 'CL', 'Location','northeast')
setfigpaper('Width',[20,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

% === DEL 
x = [1 2 3];
deltaDEL = [(OL_result.WT1.DEL_flapwise-BL_result.WT1.DEL_flapwise)/(BL_result.WT1.DEL_flapwise) (OL_result.WT2.DEL_flapwise-BL_result.WT2.DEL_flapwise)/(BL_result.WT2.DEL_flapwise) (OL_result.All.DEL_flapwise-BL_result.All.DEL_flapwise)/(BL_result.All.DEL_flapwise);
            (OL_result.WT1.DEL_edgewise-BL_result.WT1.DEL_edgewise)/(BL_result.WT1.DEL_edgewise) (OL_result.WT2.DEL_edgewise-BL_result.WT2.DEL_edgewise)/(BL_result.WT2.DEL_edgewise) (OL_result.All.DEL_edgewise-BL_result.All.DEL_edgewise)/(BL_result.All.DEL_edgewise);
            (CL_result.WT1.DEL_flapwise-BL_result.WT1.DEL_flapwise)/(BL_result.WT1.DEL_flapwise) (CL_result.WT2.DEL_flapwise-BL_result.WT2.DEL_flapwise)/(BL_result.WT2.DEL_flapwise) (CL_result.All.DEL_flapwise-BL_result.All.DEL_flapwise)/(BL_result.All.DEL_flapwise);
            (CL_result.WT1.DEL_edgewise-BL_result.WT1.DEL_edgewise)/(BL_result.WT1.DEL_edgewise) (CL_result.WT2.DEL_edgewise-BL_result.WT2.DEL_edgewise)/(BL_result.WT2.DEL_edgewise) (CL_result.All.DEL_edgewise-BL_result.All.DEL_edgewise)/(BL_result.All.DEL_edgewise)];
figure('Name', 'Power', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
bar(x,deltaDEL);
xticks(x); 
xticklabels({'WT1', 'WT2', 'WT1+WT2'}); 
ylabel('\Delta DEL [%]')
legend('OL Flapwise', 'OL Edgewise','CL Flapwise', 'CL Edgewise', 'Location','northeast')
setfigpaper('Width',[20,0.5],'Interpreter','tex','FontSize',Font,'linewidth',lw)

% === PBD
% PBD is only compared between the open-loop and closed-loop because
% baseline does not have PBD since Helix is not activated
deltaPBD = mean(OL_result.WT1.PBD) - mean(CL_result.WT1.PBD);
disp(deltaPBD)