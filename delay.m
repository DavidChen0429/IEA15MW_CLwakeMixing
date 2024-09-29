clear
close all
addpath('.\Functions');

%% Load model
original_sys = load('Model\ModelOrder4.mat');
decoupled_sys = load('Model\ModelOrder4_decoupled.mat');

%% Basic system property
eig(decoupled_sys.decouple_sys.A)
length(decoupled_sys.decouple_sys.A)
rank(ctrb(decoupled_sys.decouple_sys.A, decoupled_sys.decouple_sys.B))    
rank(obsv(decoupled_sys.decouple_sys.A, decoupled_sys.decouple_sys.C))

%% Load data
% trainData = 'train_120min_1bw_noise2.mat';       % train set
testData = 'stepResponse3.mat';                % test set
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
G = tf(decoupled_sys.decouple_sys);
delayed_sys = z^(-DeadtimeDelay) .* G;
delayed_sys = ss(delayed_sys);

%% Performance comparison 
yi2 = lsim(tf(original_sys.OLi),us2,t_test); % orignal system
yi2d = lsim(tf(decoupled_sys.decouple_sys),us2,t_test); % ss decouple 
yi2dl = lsim(delayed_sys,us2,t_test); % ss decouple 

figure()
subplot(3, 1, 1)
plot((1:length(us2)) * timeStep, us2(1, :))
hold on 
plot((1:length(us2)) * timeStep, us2(2, :))
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('\beta_{tilt}', '\beta_{yaw}')
title('Decouple Result -- Input')

subplot(3, 1, 2)
plot((1:length(yi2)) * timeStep, yi2)
hold on
plot((1:length(yi2d)) * timeStep, yi2d)
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('tilt','yaw','tilt_{dcpl}','yaw_{dcpl}')
title('Decouple Result -- Output')

subplot(3, 1, 3)
plot((1:length(yi2d)) * timeStep, yi2d)
hold on
plot((1:length(yi2dl)) * timeStep, yi2dl)
% plot(delayseq(yi2d, DeadtimeDelay))
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('tilt_{dcpl}','yaw_{dcpl}','tilt_{dly}','yaw_{dly}')
title('Decouple Result -- Output')

%% Save model
% save('Model\ModelOrder4_decoupled_delayed.mat', 'delayed_sys');