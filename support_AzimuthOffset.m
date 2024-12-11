%% Derive the Azimuth Offset of the system for decoupling
clear
close all
addpath('.\Functions');

%% Get Training data and Testing data
trainData = '120min_1bw_noise3_0AO';                    
turbineName = '.\Data\NREL5MW\';
caseName = 'sysIDE\';
IDEdata_train = load([turbineName caseName trainData]);
timeStep = 0.1;
bw = 0.0175;    % Hz
fe = 0.3*10/126;
fe = 0;
lw = 2; % line width

%% Data processing
u_train = IDEdata_train.HF_beta;
y_train = IDEdata_train.HF_helixCenter_filtered;
% Remove first few data
shiftNum = 999;
u_train = u_train(shiftNum:end, :);
y_train = y_train(shiftNum:end, :);
% Time shift the signal
DeadtimeDelay = 110; % 110 112 This should identified a semi-delayed system
u_train = u_train(1:end-DeadtimeDelay, :);
y_train = y_train(DeadtimeDelay+1:end, :);
N_train = length(u_train);
t_train = (0:N_train-1) * timeStep;
% Detrend data
u_train = detrend(u_train, 'constant');
y_train = detrend(y_train, 'constant');
% Signal scaling
us = u_train';  % 2*N
ys = y_train';  % 2*N

%% PBSID-varx System IDE
n_varx = 4;
f_varx = 200;    
p_varx = 200;

[S,X] = dordvarx(us,ys,f_varx,p_varx,'tikh','gcv');
x = dmodx(X,n_varx);

% State-Space Model Acquisition
[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,us,ys,f_varx,p_varx);
OLi = ss(Ai,Bi,Ci,Di,timeStep);
OLi.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLi.OutputName = {'z_e','y_e'};

%% Validation
% VAF
yi = lsim(OLi,us,t_train);
disp('=================================================')
disp('[Training] VAF with PBSID-varx (open loop)')
vaf(ys, yi)   

% Frequency Response
[Ga,ws] = spa_avf(us,ys,timeStep,6,[],[],'hamming');
Ga = frd(Ga,ws);
figure('Name', 'Frequence Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% bodemag(OLi, Ga);
bode(Ga, OLi);
hold on
axesHandles = findall(gcf, 'Type', 'axes');
yline(axesHandles(3), 0, 'k--', 'LineWidth', lw);
yline(axesHandles(5), 0, 'k--', 'LineWidth', lw);
yline(axesHandles(7), 0, 'k--', 'LineWidth', lw);
yline(axesHandles(9), 0, 'k--', 'LineWidth', lw);
yline(axesHandles(2), 0, 'k--', 'LineWidth', lw);
yline(axesHandles(4), 0, 'k--', 'LineWidth', lw);
yline(axesHandles(6), 0, 'k--', 'LineWidth', lw);
yline(axesHandles(8), 0, 'k--', 'LineWidth', lw);
% Don't forget to convert to Hz when using below to show bw !!!
xline(axesHandles(3), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(5), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(7), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(9), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(2), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(4), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(6), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(8), bw, 'k--', 'LineWidth', lw);
hold off;
grid on;
legend('Real','Sys-IDE','f_{bw}','Location','southeast');
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',lw)

%% Derive Azimuth Offset
G = tf(OLi);
G_resp = evalfr(G, exp(j*fe*2*pi*timeStep));
mag = abs(G_resp);
angle_rad  = angle(G_resp);
numerator = mag(1,1) * sin(angle_rad(1,1)) + mag(1,2) * cos(angle_rad(1,2));
denominator = mag(2,1) * sin(angle_rad(2,1)) + mag(2,2) * cos(angle_rad(2,2));
psi_d = numerator / denominator;
psi_0 = -psi_d;
disp(psi_0)