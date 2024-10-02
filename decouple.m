clear
close all
addpath('.\Functions');

%% Load model
buf_sys = load('Model\ModelOrder4.mat');
decouple_sys = load('Model\ModelOrder4_AzimuthOffset.mat');
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
trainData = 'train_120min_1bw_noise3%_AzimuthOffset.mat';       % train set
testData = 'stepResponse_tiltOnly_AzimuthOffset.mat';                % test set
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
u_test = detrend(u_test, 'constant');
y_test = detrend(y_test, 'constant');

% Signal scaling
us = u_train';  % 2*N
ys = y_train';  % 2*N
us2 = u_test';  % 2*N
ys2 = y_test';  % 2*N

%% Decouple 
% Original system's RGA
G_ss = dcgain(G);   % steady-state gain matrix
RGA = G_ss .* (inv(G_ss))';
eigG = eig(G_ss);

% % 1. steady-state decoupling
% ss_compensator = eye(2)/(G_ss) * 5; % calibrated for gain *5
% decouple_sys = series(ss_compensator, buf_sys.OLi);
% G_ss_dcpl = dcgain(decouple_sys);
% RGA_ssdcpl = G_ss_dcpl .* (inv(G_ss_dcpl))';

% % 2. dynamic decoupling
% % This currently has the issue of not being a non-minimal phase system
% % which introduces unstable RHP-poles when inverse G
% tzero(G)
% dyn_compensator = inv(G);
% decouple_sys2 = series(dyn_compensator, buf_sys.OLi);
% G_ss2 = dcgain(decouple_sys2);
% RGA2 = G_ss2 .* (inv(G_ss2))';

yi2 = lsim(buf_sys.OLi,us2,t_test); % testing set
yi2d = lsim(decouple_sys.OLi,us2,t_test); % ss decouple 
% yi2d2 = lsim(decouple_sys2,us2,t_test); % dyn decouple 

% yi2 = lsim(buf_sys.OLi,us,t_train); % training set
% yi2d = lsim(decouple_sys,us,t_train); % ss decouple 
% yi2d = lsim(decouple_sys2,us2,t_test); % dyn decouple 

figure()
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
plot((1:length(yi2d)) * timeStep, yi2d(:, 1),'m', 'LineWidth', 1)
hold on
plot((1:length(yi2d)) * timeStep, yi2d(:, 2),'b', 'LineWidth', 1)
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('z_e - tilt','y_e - yaw')
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
% save('Model\ModelOrder4_decoupled.mat', 'decouple_sys');