clear
close all
addpath('.\Functions');

%% Load model
buf_sys = load('Model/RightTransform_Azimuth96/ModelOrder4_noise1p_opposite.mat');
% decouple_sys = load('Model\ModelOrder4_AzimuthOffset.mat');
A = buf_sys.OLi.A;
B = buf_sys.OLi.B;
C = buf_sys.OLi.C;
D = buf_sys.OLi.D;

sys = buf_sys.OLi;
G = tf(buf_sys.OLi);        % transfer matrix 

%% Basic system property
fprintf('======== System property \n');
fprintf(' System Dimension: %.0f \n', size(A, 1));
fprintf(' Ctrb Matrix Rank: %.0f \n', rank(ctrb(A, B)));
fprintf(' Obsv Matrix Rank: %.0f \n', rank(obsv(A, C)));
fprintf(' Eigenvalues of A: \n');
disp(eig(A));

%% Load data
trainData = 'train_240min_1bw_noise1%_AzimuthOffset96.mat'; 
testData = 'stepResponse_yawOnly_AzimuthOffset96.mat';    
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\sysIDE\';
IDEdata_train = load([turbineName caseName trainData]);
IDEdata_test = load([turbineName caseName testData]);
timeStep = 0.1;

u_train = IDEdata_train.HF_beta;
y_train = IDEdata_train.HF_helixCenter_filtered;
u_test = IDEdata_test.HF_beta;
y_test = IDEdata_test.HF_helixCenter_filtered;

% Remove first few data
shiftNum = 500;
u_train = u_train(shiftNum:end, :);
y_train = y_train(shiftNum:end, :);
u_test = u_test(shiftNum:end, :);
y_test = y_test(shiftNum:end, :);

% Time shift the signal
DeadtimeDelay = 110;
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
% u_test = detrend(u_test, 'constant');
% y_test = detrend(y_test, 'constant');

% Signal scaling
us = u_train';  % 2*N
ys = y_train';  % 2*N
us2 = u_test';  % 2*N
ys2 = y_test';  % 2*N

% for the oppposite system
ys(1, :) = -1 * ys(1, :);
ys2(1, :) = -1 * ys2(1, :);

%% Decouple 
% Original system's RGA
G_ss = dcgain(G);   % steady-state gain matrix
RGA = G_ss .* (inv(G_ss))';
eigG = eig(G_ss);
disp(RGA)
bwG11 = calculateBandwidth(G(1, 1));
bwG22 = calculateBandwidth(G(2, 2));
bw = (bwG11 + bwG22)/2; 
bw = 0.0175;
testCoupling(buf_sys.OLi, bw, timeStep);

% % 1. steady-state decoupling
ss_compensator = eye(2)/(G_ss) * 3; % match gain *3
OLi = series(ss_compensator, buf_sys.OLi);
OLi.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLi.OutputName = {'z_e','y_e'};
testCoupling(OLi, bw, timeStep);

% 2. Specific frequency decouple / pre-compensator
G_bw = evalfr(sys, exp(j*bw*2*pi*timeStep));
G_bw_real = abs(G_bw);
bw_compensator = inv(G_bw_real);  % * 3 
OLi2 = series(bw_compensator, buf_sys.OLi);
OLi2.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLi2.OutputName = {'z_e','y_e'};
testCoupling(OLi2, bw, timeStep);

% 3. dynamic decoupling
% This currently has the issue of not being a non-minimal phase system
% which introduces unstable RHP-poles when inverse G, break the stability
% tzero(G)
dyn_compensator = inv(G);
OLi3 = series(dyn_compensator, buf_sys.OLi);
OLi3.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLi3.OutputName = {'z_e','y_e'};
testCoupling(OLi3, bw, timeStep);

%% Compare result
yi2_train = lsim(buf_sys.OLi,us,t_train); % testing set
% yi2d_train = lsim(OLi,us,t_train); % ss decouple 
yi2d_train = lsim(OLi2,us,t_train); % bw decouple 
yi2_test = lsim(buf_sys.OLi,us2,t_test); % testing set
% yi2d_test = lsim(OLi,us2,t_test); % ss decouple 
yi2d_test = lsim(OLi2,us2,t_test); % bw decouple 

% % VAF
% disp('=================================================')
% disp('[Training] VAF with PBSID-varx (open loop)')
% vaf(ys, yi2_train)
% vaf(ys, yi2d_train)
% vaf(ys, yi2d2_train)
% disp('[Testing] VAF with PBSID-varx (open loop)')
% vaf(ys2, yi2_test)  
% vaf(ys2, yi2d_test) 
% vaf(ys2, yi2d2_test)  

%% Frequency Domain Fitting Result
lw = 1;
[Ga,ws] = spa_avf(us,ys,timeStep,6,[],[],'hamming');
Ga = frd(Ga,ws);
figure('Name', 'Frequence Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% bodemag(buf_sys.OLi, OLi, OLi2, Ga);
% bode(Ga, buf_sys.OLi, OLi, OLi2);
bode(buf_sys.OLi, OLi, OLi2);
hold on
axesHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(axesHandles)
    yline(axesHandles(k), 0, 'k--', 'LineWidth', lw);  
end
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
grid on
% legend('Real', 'Original Sys','Ss', 'Bw','Location','southeast');
legend('Original','Steady-state', 'Bandwidth','Location','southeast');
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',15,'linewidth',lw)

[U,S,V] = svd(abs(G_bw));

%% Test Set
lw = 1;
figure('Name', 'Time-domain Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot((1:length(us2)) * timeStep, us2(1, :), 'm', 'LineWidth', lw)
hold on 
plot((1:length(us2)) * timeStep, us2(2, :), 'b', 'LineWidth', lw)
yline(0, '--', 'LineWidth', lw)
hold off
xlabel('Time [s]')
xlim([0, length(us2)*timeStep])
ylabel('Magnitude [m]')
legend('\beta^e_{tilt}', '\beta^e_{yaw}','Location','southeast')
title('Decouple Result -- Input')

% For already decoupled system
subplot(2, 1, 2)
plot((1:length(yi2_test)) * timeStep, yi2_test(:, 1),'m--', 'LineWidth', lw)
hold on
plot((1:length(yi2_test)) * timeStep, yi2_test(:, 2),'b--', 'LineWidth', lw)
plot((1:length(yi2d_test)) * timeStep, yi2d_test(:, 1),'m', 'LineWidth', lw)
plot((1:length(yi2d_test)) * timeStep, yi2d_test(:, 2),'b', 'LineWidth', lw)
yline(0, '--', 'LineWidth', lw)
hold off
xlabel('Time [s]')
xlim([0, length(yi2_test)*timeStep])
ylabel('Magnitude [m]')
legend('z^e','y^e','z^e_d','y^e_d','Location','southeast')
title('Decouple Result -- Output')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',lw)