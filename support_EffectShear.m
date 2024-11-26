clear
% close all 
addpath('.\Functions');
%clc

%% Data file (Chage this accordingly)
turbineName = '.\Data\NREL5MW\';
caseName = 'Experiment\Str0.3_U10_1Dd_10Hz_CCW\2TurbinesNew\';

% Different case
windCase = 'Shear2'; % Shear, Turb, Both, Shear2
basefileOL = '2Turbines_OL_Helix_mag3_4D.mat';

if strcmp(windCase, 'Shear')
    Krakenfile = '2Turbines_Baseline_Shear0.2_4D.mat';
    KrakenOLfile = '2Turbines_OL_Helix_Shear0.2_mag3_4D.mat';
elseif strcmp(windCase, 'Shear2')
    Krakenfile = '2Turbines_Baseline_Shear0.3_4D.mat';
    KrakenOLfile = '2Turbines_OL_Helix_Shear0.3_mag3_4D.mat';
end

BasicOL = load([turbineName caseName basefileOL]);
Kraken = load([turbineName caseName Krakenfile]);
KrakenOL = load([turbineName caseName KrakenOLfile]);

simLength = length(BasicOL.Power_store);
timeStep = 0.1;
filter = 3000;
t = (1:(simLength-filter+1)) * timeStep;

%% Visualization
% Hub jet in shear (HF)
figure('Name', 'Hub Jet in Shear HF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t, BasicOL.HF_helixCenter_filtered(filter:end, 1))
hold on
plot(t, KrakenOL.HF_helixCenter_filtered(filter:end, 1))
hold off
title('z^e')
ylim([0 10])
legend('Uniform', 'Shear')
subplot(2, 1, 2)
plot(t, BasicOL.HF_helixCenter_filtered(filter:end, 2))
hold on
plot(t, KrakenOL.HF_helixCenter_filtered(filter:end, 2))
hold off
title('y^e')
ylim([0 10])
legend('Uniform', 'Shear')

% Hub jet in shear (FF)
figure('Name', 'Hub Jet in Shear HF', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t, BasicOL.FF_helixCenter_filtered(filter:end, 1))
hold on
plot(t, KrakenOL.FF_helixCenter_filtered(filter:end, 1))
hold off
title('z^e')
legend('Uniform', 'Shear')
subplot(2, 1, 2)
plot(t, BasicOL.FF_helixCenter_filtered(filter:end, 2))
hold on
plot(t, KrakenOL.FF_helixCenter_filtered(filter:end, 2))
hold off
title('y^e')
legend('Uniform', 'Shear')