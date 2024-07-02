%% Helix of IEA15MW in script
clear
close all 
%clc

%% Define paths
UserPath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\IEA15MW_CLwakeMixing\'; 
% QBladePath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\QBladeCE_2.0.7.1\'; 
QBladePath = 'C:\Users\DAVID CHEN\Desktop\TU_Delft\Thesis\QBladeEE_2.0.6.4\'; 
SourcePath = [UserPath 'Source\'];
% DllPath = [QBladePath 'QBladeCE_2.0.7.dll'];
DllPath = [QBladePath 'QBladeEE_2.0.6.dll'];
simFile = [SourcePath 'IEA15MW_torque_Helix.sim'];
addpath('.\Functions');

% loadlibrary(DllPath,'.\QBladeDLLFunctions.h','alias','QBladeDLL') 
loadlibrary(DllPath,'QBladeLibInclude.h','alias','QBladeDLL') 
m = libfunctions('QBladeDLL');
if isempty(m)
    fprintf('Error')
end

% unloadlibrary 'QBladeDLL'