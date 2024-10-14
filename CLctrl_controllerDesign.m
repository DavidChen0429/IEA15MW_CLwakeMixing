% This file contains controller design for the decoupled MIMO
% system based on the identified decoupled model sys. This work
% based on the observation that y plays a dominate role in overall ouput

clear
close all
addpath('.\Functions');

%% Load model
buf_sys = load('Model\RightTransform_Azimuth96\ModelOrder4_noise1p_opposite_decoupled.mat');
A = buf_sys.OLi.A;
B = buf_sys.OLi.B;
C = buf_sys.OLi.C;
D = buf_sys.OLi.D;

sys = buf_sys.OLi;
sys.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
sys.OutputName = {'z_e','y_e'};
% G = tf(buf_sys.OLi);        % transfer matrix 
G = tf(buf_sys.OLi);        % transfer matrix 
% G = G([2, 1], :);

% Basic system property
eig(A)
size(A, 1)
rank(ctrb(A, B))    
rank(obsv(A, C))

% Load data
trainData = 'train_120min_1bw_noise5%_AzimuthOffset6.mat';       % train set
testData = 'stepResponse_both_AzimuthOffset6.mat';                % test set
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
shiftNum = 999;
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
% [us,Du,ys,Dy] = sigscale(u_train, y_train);
% [us2,Du2,ys2,Dy2] = sigscale(u_test, y_test);
us = u_train';  % 2*N
ys = y_train';  % 2*N
us2 = u_test';  % 2*N
ys2 = y_test';  % 2*N

% Bode diagram (Frequency domain response)
figure('Name', 'Bode Diagram OL System', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
bode(G);
title('Bode Diagram of Open-Loop System')

%% PID Controller Design
% func: pidtune
wc = 0.010;
C11 = pidtune(G(1,1), 'PI', wc);
C22 = pidtune(G(2,2), 'PI', wc);
C12 = 0;
C21 = 0;
% Kp = 0; % 0
% Ki = 0.0375; % 0.05
% Ts = timeStep;
% C11 = pid(Kp, Ki, 0, 0, Ts);
C_mimo = [C11, 0;
          0, 0];
OL_ctrl = C_mimo * G;

%%%
% Open loop transfer matrix becomes
%   [C11*G11 C11*G12]
%   [C22*G21 C22*G22]
%%%

% % Before control
% % C11 in charge of below
% figure('Name', 'G11 G12 BeforeCtrl', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% subplot(1,2,1)
% margin(G(1, 1));
% subplot(1,2,2)
% margin(G(1, 2)); 
% % C22 in charge of below
% figure('Name', 'G21 G22 BeforeCtrl', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% subplot(1,2,1)
% margin(G(2, 1)); 
% subplot(1,2,2)
% margin(G(2, 2));
% 
% % Loop shaping
% % bode(OL_ctrl);
% % C11 in charge of below
% figure('Name', 'G11 G12 AfterCtrl', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% subplot(1,2,1)
% margin(OL_ctrl(1, 1));
% subplot(1,2,2)
% margin(OL_ctrl(1, 2)); 
% % C22 in charge of below
% figure('Name', 'G21 G22 AfterCtrl', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% subplot(1,2,1) 
% margin(OL_ctrl(2, 1)); 
% subplot(1,2,2)
% margin(OL_ctrl(2, 2));

closed_loop_sys = feedback(OL_ctrl, eye(2));
t = 0:timeStep:1000;  % Time vector for simulation
figure('Name', 'After Control CL Step', 'NumberTitle', 'off');
step(closed_loop_sys, t);
title('Controlled CL System');
grid on;

%% Faster tuning
% close all
Kp = 0; % 0
Ki = 0.0375; % 0.0375
Ts = timeStep;
C22 = pidtune(G(2,2), 'I', 0.010);
C11 = pid(Kp, Ki, 0, 0, Ts);
C_mimo = [C11, 0;
          0, C22];
OL_ctrl = C_mimo * G;

% C11 in charge of below
figure('Name', 'G11 G12 AfterCtrl', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(1,2,1) 
margin(OL_ctrl(1, 1)); 
subplot(1,2,2)
margin(OL_ctrl(1, 2)); % Fucked

% % C22 in charge of below
% figure('Name', 'G21 G22 AfterCtrl', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% subplot(1,2,1) 
% margin(OL_ctrl(2, 1)); % Fucked
% subplot(1,2,2)
% margin(OL_ctrl(2, 2));

%% Phase-lead compensator design
zpk(G(1,1))
zpk(G(1,2))
zpk(G(2,1))
zpk(G(2,2))