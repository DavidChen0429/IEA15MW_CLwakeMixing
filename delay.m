clear
close all
addpath('.\Functions');

%% Load model
original_sys = load('Model\ModelOrder4.mat');
decoupled_sys = load('Model\ModelOrder4_AzimuthOffset.mat');

%% Basic system property
eig(decoupled_sys.OLi.A)
length(decoupled_sys.OLi.A)
rank(ctrb(decoupled_sys.OLi.A, decoupled_sys.OLi.B))    
rank(obsv(decoupled_sys.OLi.A, decoupled_sys.OLi.C))

%% Load data
% trainData = 'train_120min_1bw_noise2.mat';       % train set
testData = 'stepResponse_both_AzimuthOffset.mat';                % test set
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\sysIDE\';
% IDEdata_train = load([turbineName caseName trainData]);
IDEdata_test = load([turbineName caseName testData]);
timeStep = 0.1;

% u_train = IDEdata_train.HF_beta;
% y_train = IDEdata_train.HF_helixCenter_filtered;
u_test = IDEdata_test.HF_beta;
y_test = IDEdata_test.HF_helixCenter_filtered;

% Remove first few data
shiftNum = 800;
% u_train = u_train(shiftNum:end, :);
% y_train = y_train(shiftNum:end, :);
u_test = u_test(shiftNum:end, :);
y_test = y_test(shiftNum:end, :);

% Time shift the signal
DeadtimeDelay = 110;
% u_train = u_train(1:end-DeadtimeDelay, :);
% y_train = y_train(DeadtimeDelay+1:end, :);
% N_train = length(u_train);
% t_train = (0:N_train-1) * timeStep;

u_test = u_test(1:end-DeadtimeDelay, :);
y_test = y_test(DeadtimeDelay+1:end, :);
N_test = length(u_test);
t_test = (0:N_test-1) * timeStep;

% Detrend data
% u_train = detrend(u_train, 'constant');
% y_train = detrend(y_train, 'constant');
u_test = detrend(u_test, 'constant');
y_test = detrend(y_test, 'constant');

% Signal scaling
% [us,Du,ys,Dy] = sigscale(u_train, y_train);
% [us2,Du2,ys2,Dy2] = sigscale(u_test, y_test);
% us = u_train';  % 2*N
% ys = y_train';  % 2*N
us2 = u_test';  % 2*N
ys2 = y_test';  % 2*N

%% Add delay factor
% DeadtimeDelay = 1;
% A = decoupled_sys.decouple_sys.A;
% B = decoupled_sys.decouple_sys.B;
% C = decoupled_sys.decouple_sys.C;
% D = decoupled_sys.decouple_sys.D;
% A_aug = [A, zeros(size(A,1), size(A,1)*DeadtimeDelay);
%          zeros(size(A,1)*DeadtimeDelay, size(A,1)), eye(size(A,1)*DeadtimeDelay)];
% B_aug = [zeros(size(A, 1)*DeadtimeDelay, 2); 
%          B];
% C_aug = [zeros(size(C,1)*DeadtimeDelay, size(C,2)), eye(size(C,1)*DeadtimeDelay);
%          C, zeros(size(C,1), size(C,2)*DeadtimeDelay)];
% D_aug = D;
% delayed_sys_ss = ss(A_aug, B_aug, C_aug, D_aug);

z = tf('z', timeStep);
G = tf(decoupled_sys.OLi);
OLi = z^(-DeadtimeDelay) .* G;
OLi = ss(OLi);

G2 = tf(original_sys.OLi);
delayed_sys_original = z^(-DeadtimeDelay) .* G2;
delayed_sys_original = ss(delayed_sys_original);

%% Performance comparison 
yi2 = lsim(tf(original_sys.OLi),us2,t_test); % orignal system
yi2d = lsim(tf(decoupled_sys.OLi),us2,t_test); % ss decouple 
yi2dl = lsim(OLi,us2,t_test); % ss decouple 

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
title('Input')

subplot(2, 1, 2)
plot((1:length(yi2d)) * timeStep, yi2d(:, 1), 'm', 'LineWidth', 1)
hold on
plot((1:length(yi2d)) * timeStep, yi2d(:, 2), 'b', 'LineWidth', 1)
plot((1:length(yi2dl)) * timeStep, yi2dl(:, 1), 'm--', 'LineWidth', 1)
plot((1:length(yi2dl)) * timeStep, yi2dl(:, 2), 'b--', 'LineWidth', 1)
% plot(delayseq(yi2d, DeadtimeDelay))
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('z_e','y_e','z_{e,d}','y_{e,d}')
title('Delayed Result -- Output')

%% Save model
% save('Model\ModelOrder4_AzimuthOffset_delayed.mat', 'OLi');
% save('Model\ModelOrder4_delayed.mat', 'delayed_sys_original');