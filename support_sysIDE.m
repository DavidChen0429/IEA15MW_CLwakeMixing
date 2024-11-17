clear
close all
addpath('.\Functions');

%% Get Training data and Testing data
trainData = 'train_240min_1bw_noise1%_AzimuthOffset96.mat';      
testData = 'stepResponse_both_AzimuthOffset96.mat';                
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\sysIDE\';
IDEdata_train = load([turbineName caseName trainData]);
IDEdata_test = load([turbineName caseName testData]);
timeStep = 0.1;
bw = 0.0175;    % Hz
lw = 2; % line width

u_train = IDEdata_train.HF_beta;
y_train = IDEdata_train.HF_helixCenter_filtered;
u_test = IDEdata_test.HF_beta;
y_test = IDEdata_test.HF_helixCenter_filtered;

% Remove first few data
shiftNum = 999;
u_train = u_train(shiftNum:end, :);
y_train = y_train(shiftNum:end, :);
u_test = u_test(shiftNum:end, :);
y_test = y_test(shiftNum:end, :);

% Time shift the signal
DeadtimeDelay = 110; % 110 112 This should identified a semi-delayed system
u_train = u_train(1:end-DeadtimeDelay, :);
y_train = y_train(DeadtimeDelay+1:end, :);
N_train = length(u_train);
t_train = (0:N_train-1) * timeStep;

u_test = u_test(1:end-DeadtimeDelay, :);
y_test = y_test(DeadtimeDelay+1:end, :);
N_test = length(u_test);
t_test = (0:N_test-1) * timeStep;

% Detrend data
u_train = detrend(u_train, 'constant');
y_train = detrend(y_train, 'constant');
u_test = detrend(u_test, 'constant');
y_test = detrend(y_test, 'constant');

% Signal scaling
% [us,Du,ys,Dy] = sigscale(u_train, y_train);
% [us2,Du2,ys2,Dy2] = sigscale(u_test, y_test);
us = u_train';  % 2*N
ys = y_train';  % 2*N
us2 = u_test';  % 2*N
ys2 = y_test';  % 2*N

ys(1, :) = -1 * ys(1, :);
ys2(1, :) = -1 * ys2(1, :);

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

% Validation (VAF)
yi = lsim(OLi,us,t_train);
disp('=================================================')
disp('[Training] VAF with PBSID-varx (open loop)')
vaf(ys, yi)   
yi_test = lsim(OLi,us2,t_test);
disp('[Testing] VAF with PBSID-varx (open loop)')
vaf(ys2, yi_test)  

%% Power Spectrum Density (Input Signal)
Ts_prbn = timeStep;
[M1,F1] = pwelch(us(1, :),[],[],[],1/Ts_prbn);
[M2,F2] = pwelch(us(2, :),[],[],[],1/Ts_prbn);
figure('Name', 'Input PSD', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(1,2,1)
semilogx(F1,mag2db(M1),'k','LineWidth',lw)
hold on
semilogx(F2,mag2db(M2),'r','LineWidth',lw)
yline(0, ':', 'LineWidth', lw)
xline(bw, '--', 'LineWidth', lw)
hold off
grid on
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');
legend('\beta^e_{tilt}', '\beta^e_{yaw}','f_{bw}', 'Location', 'southeast')
title('Input PSD')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',lw)

[M1,F1] = pwelch(ys(1, :),[],[],[],1/Ts_prbn);
[M2,F2] = pwelch(ys(2, :),[],[],[],1/Ts_prbn);
% figure('Name', 'Output PSD', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(1,2,2)
semilogx(F1,mag2db(M1),'k','LineWidth',lw)
hold on
grid on
semilogx(F2,mag2db(M2),'r','LineWidth',lw)
yline(0, ':', 'LineWidth', lw)
xline(bw, '--', 'LineWidth', lw)
hold off
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');
legend('z_f', 'y_f','f_{bw}', 'Location', 'southeast')
title('Output PSD')

setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',lw)

%% Singular Value Visualization
figure('Name', 'Singular Value Visualization', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
semilogy(S,'*');
title('Singular Value PBSID-opt')
xlim([0 50]);
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',lw)

%% Frequency Domain Fitting Result
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

%% Time Domain Fitting Result (Only testset)
yi2 = lsim(OLi,us2,t_test); % testing set
figure('Name', 'Time Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2,1,1)
plot((1:length(yi2)) * timeStep, yi2(:, 1), 'LineWidth', lw)
hold on
plot((1:length(yi2)) * timeStep, ys2(1, :), 'LineWidth', lw)
plot((1:length(us2)) * timeStep, us2(1, :), 'k:', 'LineWidth', lw)
yline(0, '--', 'LineWidth', 1)
hold off
legend('Sys-IDE', 'Real','\beta^e_{tilt}', 'Location', 'southeast')
xlabel('Time [s]')
ylabel('z_e [m]')
title('Test Set')
subplot(2,1,2)
plot((1:length(yi2)) * timeStep, yi2(:, 2), 'LineWidth', lw)
hold on
plot((1:length(yi2)) * timeStep, ys2(2, :), 'LineWidth', lw)
plot((1:length(us2)) * timeStep, us2(2, :), 'k:', 'LineWidth', lw)
yline(0, '--', 'LineWidth', 1)
hold off
legend('Sys-IDE', 'Real','\beta^e_{yaw}', 'Location', 'southeast')
xlabel('Time [s]')
ylabel('y_e [m]')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',lw)

%% Numerical Result Check
% ===== Coupling
% ss
G = tf(OLi);
G_ss = dcgain(G);   % steady-state gain matrix
RGA = G_ss .* (inv(G_ss))';
disp(RGA);
% bw
G_bw = evalfr(G, exp(j*bw*2*pi*timeStep));
G_bw_real = abs(G_bw);
RGA2 = G_bw_real .* (inv(G_bw_real))';
disp(RGA2);

% bandwidth
calculateBandwidth(G(1,1));
calculateBandwidth(G(2,2));