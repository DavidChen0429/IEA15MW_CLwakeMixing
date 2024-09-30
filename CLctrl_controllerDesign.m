% This file contains controller design for the decoupled MIMO
% system based on the identified decoupled model sys. This work
% based on the observation that y plays a dominate role in overall ouput

clear
close all
addpath('.\Functions');

%% Load model
buf_sys = load('Model\ModelOrder4_decoupled.mat');
A = buf_sys.decouple_sys.A;
B = buf_sys.decouple_sys.B;
C = buf_sys.decouple_sys.C;
D = buf_sys.decouple_sys.D;

sys = buf_sys.decouple_sys;
sys.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
sys.OutputName = {'z_e','y_e'};
G = tf(buf_sys.decouple_sys);        % transfer matrix 

%% Basic system property
eig(A)
length(A)
rank(ctrb(A, B))    
rank(obsv(A, C))

%% Load data
trainData = 'train_120min_1bw_noise2.mat';       % train set
testData = 'stepResponse3.mat';                % test set
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
% [us,Du,ys,Dy] = sigscale(u_train, y_train);
% [us2,Du2,ys2,Dy2] = sigscale(u_test, y_test);
us = u_train';  % 2*N
ys = y_train';  % 2*N
us2 = u_test';  % 2*N
ys2 = y_test';  % 2*N

%% Step response (Time domain response)
% open loop (stable)
figure;
step(sys);
title('Step Response of the OL System');
xlabel('Time (seconds)');
ylabel('Response');

% closed loop (unstable)
figure
CL_noControl = feedback(sys, eye(2));
step(feedback(sys, eye(2)))
title('Step Response of the CL System');
xlabel('Time (seconds)');
ylabel('Response');

%% Bode diagram (Frequency domain response)
figure;
bode(sys);
title('Bode Diagram of Open-Loop System')

figure
margin(G(1, 1));
figure
margin(G(1, 2));
figure
margin(G(2, 1));
figure
margin(G(2, 2));

%% PID Controller Design
% func: pidtune
% C11 = pidtune(G(1,1), 'PI');
% C22 = pidtune(G(2,2), 'PI');
C11 = pid(1,0,0,0,timeStep);
C22 = pid(1,0,0,0,timeStep);
C_mimo = [C11, 0;
          0, C22];

closed_loop_sys = feedback(C_mimo * G, eye(2));
t = 0:timeStep:10;  % Time vector for simulation
step(closed_loop_sys, t);
title('Closed-Loop Response with Diagonal MIMO PID Controller');
grid on;