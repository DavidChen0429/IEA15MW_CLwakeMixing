%%
clear
close all 
clc

% Define paths
UserPath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\IEA15MW_Helix_LiDAR\'; 
QBladePath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\QBladeEE_2.0.6.4\'; 
disconPath = QBladePath;
MatlabPath = UserPath;
SourcePath = [UserPath 'Source'];
DllPath = [QBladePath 'QBladeEE_2.0.6.dll'];
addpath('.\Functions');

loadlibrary(DllPath,'.\QBladeDLLFunctions.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');

if isempty(m)
    fprintf('Error')
end

% Defining Settings
Hub_NREL5MW = 90;
Str = 0.4;
D_NREL5MW = 126;                            % rotor diameter NREL5MW
U_NREL5MW = 9;                      % inflow velocity NREL5MW
Helix_amplitude = 4;                

IPC = 0;
Pulse = 0/57.3;                     % 1 rad = 57.3 degree
Helix_CM2 = Helix_amplitude/57.3;
Helix_CM3 = Helix_amplitude/57.3;
Freq = Str*U_NREL5MW/D_NREL5MW*2*pi ;       % Definition of Strouhal number  
Pitch_off = 0/57.3;

SimNr = length(Freq);
tic

%%
for i_for1 = 1:1:SimNr

    Sim_name_folder = ['Sim_', num2str(i_for1)];
    
    copyfile([SourcePath],Sim_name_folder);

    cd(disconPath)

    writeDisconIPC(IPC,Pulse,Helix_CM2,Helix_CM3,Freq(i_for1),Pitch_off)
    cd(MatlabPath)
    
    run IPC_Simulation.m
    
    cd(Sim_name_folder)
    
    save(['Output_',num2str(i_for1)],'PitchAngles')
    %save(['Output_',num2str(i_for1)],'windspeed')
    
    cd(disconPath)
    delete('discon.in');
    cd(MatlabPath)
end
%%
toc

figure;
plot(PitchAngles(:,1))
hold on
plot(PitchAngles(:,2))
plot(PitchAngles(:,3))
grid on
legend('Blade 1','Blade 2','Blade 3')

%%
%run Animation.m