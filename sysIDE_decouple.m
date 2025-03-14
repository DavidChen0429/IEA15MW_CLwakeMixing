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
eig(A)
length(A)
rank(ctrb(A, B))    
rank(obsv(A, C))

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
shiftNum = 800;
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
yi2d_train = lsim(OLi,us,t_train); % ss decouple 
yi2d2_train = lsim(OLi2,us,t_train); % bw decouple 
yi2_test = lsim(buf_sys.OLi,us2,t_test); % testing set
yi2d_test = lsim(OLi,us2,t_test); % ss decouple 
yi2d2_test = lsim(OLi2,us2,t_test); % bw decouple 

% VAF
disp('=================================================')
disp('[Training] VAF with PBSID-varx (open loop)')
vaf(ys, yi2_train)
vaf(ys, yi2d_train)
vaf(ys, yi2d2_train)
disp('[Testing] VAF with PBSID-varx (open loop)')
vaf(ys2, yi2_test)  
vaf(ys2, yi2d_test) 
vaf(ys2, yi2d2_test)  

% Validation (Frequency domain)
[Ga,ws] = spa_avf(us,ys,timeStep,6,[],[],'hamming');
Ga = frd(Ga,ws);
figure('Name', 'Frequence Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
bodemag(buf_sys.OLi, OLi, OLi2, Ga, 'c:');
hold on
axesHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(axesHandles)
    yline(axesHandles(k), 0, 'k--', 'LineWidth', 0.5);  
end
xline(axesHandles(3), bw*2*pi, 'k--', 'LineWidth', 0.5);
xline(axesHandles(5), bw*2*pi, 'k--', 'LineWidth', 0.5);
xline(axesHandles(7), bw*2*pi, 'k--', 'LineWidth', 0.5);
xline(axesHandles(9), bw*2*pi, 'k--', 'LineWidth', 0.5);
% when switching to Hz as unit, use the below one
% xline(axesHandles(3), bw, 'k--', 'LineWidth', 0.5);
% xline(axesHandles(5), bw, 'k--', 'LineWidth', 0.5);
% xline(axesHandles(7), bw, 'k--', 'LineWidth', 0.5);
% xline(axesHandles(9), bw, 'k--', 'LineWidth', 0.5);
hold off;
grid on
legend('Original Sys','SS Decoupled Sys', 'BW Decoupled Sys', 'Real');

%% Test Set
% % Frequency response
% figure('Name', 'Original System', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% bode(G)
% figure('Name', 'Decoupled System', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% bode(OLi)

figure('Name', 'Time-domain Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot((1:length(us2)) * timeStep, us2(1, :), 'm', 'LineWidth', 1)
hold on 
plot((1:length(us2)) * timeStep, us2(2, :), 'b', 'LineWidth', 1)
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('\beta^e_{tilt}', '\beta^e_{yaw}')
title('Decouple Result -- Input')

% For already decoupled system
subplot(2, 1, 2)
plot((1:length(yi2_test)) * timeStep, yi2_test(:, 1),'m--', 'LineWidth', 1)
hold on
plot((1:length(yi2_test)) * timeStep, yi2_test(:, 2),'b--', 'LineWidth', 1)
plot((1:length(yi2d_test)) * timeStep, yi2d_test(:, 1),'m', 'LineWidth', 1)
plot((1:length(yi2d_test)) * timeStep, yi2d_test(:, 2),'b', 'LineWidth', 1)
plot((1:length(yi2d_test)) * timeStep, yi2d2_test(:, 1),'m:', 'LineWidth', 1)
plot((1:length(yi2d_test)) * timeStep, yi2d2_test(:, 2),'b:', 'LineWidth', 1)
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
% legend('z_e - tilt','y_e - yaw','z_{e,d} - tilt','y_{e,d} - yaw')
legend('z_e - tilt','y_e - yaw','z_{e,d} - tilt','y_{e,d} - yaw','z_{e,d2} - tilt','y_{e,d2} - yaw')
% legend('tilt','yaw','tilt_{dcpl}','yaw_{dcpl}','tilt_{dcpl2}','yaw_{dcpl2}')
title('Decouple Result -- Output')

% When you actually decouple
% subplot(2, 1, 2)
% plot((1:length(yi2)) * timeStep, yi2(:, 1),'m')
% hold on
% plot((1:length(yi2)) * timeStep, yi2(:, 2),'b')
% plot((1:length(yi2d)) * timeStep, yi2d(:, 1),'r--', 'LineWidth', 1)
% plot((1:length(yi2d)) * timeStep, yi2d(:, 2),'b--', 'LineWidth', 1)
% % plot(yi2d2)
% yline(0, '--', 'LineWidth', 1)
% hold off
% xlabel('Time [s]')
% ylabel('Magnitude')
% legend('tilt','yaw','tilt_{dcpl}','yaw_{dcpl}')
% % legend('tilt','yaw','tilt_{dcpl}','yaw_{dcpl}','tilt_{dcpl2}','yaw_{dcpl2}')
% title('Decouple Result -- Output')

% ss_compensator = [-0.7039   -0.4862;
%                   -0.4806    0.6978];

% % Question: can I really use a compensator like that ????
% figure
% buf_result = yi2 * ss_compensator;
% subplot(2, 1, 1)
% plot((1:length(yi2d)) * 0.1, yi2d(:, 1))
% hold on
% plot((1:length(buf_result)) * 0.1, buf_result(:, 1))
% yline(0, '--', 'LineWidth', 1)
% hold off
% legend('z_{d1}','z_{d2}')
% xlabel('Time [s]')
% ylabel('Magnitude')
% title("Output Comparison")
% subplot(2, 1, 2)
% plot((1:length(yi2d)) * 0.1, yi2d(:, 2))
% hold on
% plot((1:length(buf_result)) * 0.1, buf_result(:, 2))
% yline(0, '--', 'LineWidth', 1)
% hold off
% legend('y_{d1}','y_{d2}')
% xlabel('Time [s]')
% ylabel('Magnitude')

%% Frequency responses
% [Ga,ws] = spa_avf(us,ys,timeStep,6,[],[],'hamming');
% Ga = frd(Ga,ws);
% figure 
% bodemag(Ga,buf_sys.OLi, decouple_sys);
% legend('SPA','Identified','Decoupled')

%% Save model
% save('Model/RightTransform_Azimuth96/ModelOrder4_noise1p_opposite_decoupled2.mat', 'OLi');
% save('Model/RightTransform_Azimuth96/ModelOrder4_noise1p_opposite_Dyndecoupled.mat', 'OLi3');