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
bw = calculateBandwidth(G(1, 1));   % Hz
grid on
% hold on
% ax = findall(gcf, 'Type', 'axes');  % Find all axes in the current figure
% for i = 1:length(ax)
%     xline(ax(i), 2*pi*bw, 'r:', 'LineWidth', 1);
% end
% hold off
title('Bode Diagram of Open-Loop System')

%% MIMO Controller Design (H infinity)
W_s = tf([1, 1.6], [100, 1]);  % Emphasizes performance and disturbance rejection
W_t = tf([1, 1], [1, 1]);  % Emphasizes robustness and noise rejection
% Working weight function
% W_s = tf([1, 1.6], [100, 1]);
% W_t = tf([1, 1], [1, 1]);
W_s_d = c2d(W_s, timeStep, 'tustin');
W_t_d = c2d(W_t, timeStep, 'tustin');
P = augw(sys, W_s_d, [], W_t_d);  % augw creates the weighted augmented plant
ncont = 2; 
nmeas = 2; 
[K,CL,gamma] = hinfsyn(P,nmeas,ncont);
% [K,CL,gamma] = hinfsyn(sys,nmeas,ncont); % for debug code
sys_cl = feedback(sys, K);

% Step response check
figure('Name', 'Step Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
step(sys)
hold on
step(sys_cl)
hold off
legend('OL System', 'CL System')
grid on

%% Default vs. Iterative
close all
% 1. Default
closed_loop_sys = feedback(sys, K);
t = 0:timeStep:200;  % Time vector for simulation
figure('Name', 'Default Way', 'NumberTitle', 'off', 'Position', [100, 450, 400, 300]);
[y, tOut] = step(closed_loop_sys, t);
title('Controlled CL System');
plot(tOut, y(:, 1) + y(:, 3))
% plot(tOut, y(:, 1))
hold on
plot(tOut, y(:, 2) + y(:, 4))
% plot(tOut, y(:, 2))
% plot(tOut, y(:, 3))
% plot(tOut, y(:, 4))
hold off
xlabel('Time (s)');
ylabel('System Output');
title('Controlled CL System');
legend('z_e', 'y_e');
% legend('\beta^e_{tilt} - z_e', '\beta^e_{tilt} - y_e', '\beta^e_{yaw} - z_e', '\beta^e_{yaw} - y_e');
grid on;

% Iterative debug (Below code is right)
A_cl = closed_loop_sys.A;
B_cl = closed_loop_sys.B;
C_cl = closed_loop_sys.C;
D_cl = closed_loop_sys.D;
t = 0:timeStep:200;
x = zeros(length(A_cl), length(t)+1);
y = zeros(length(D_cl), length(t)); 
r = 1*ones(length(B(1, :)), length(t)+1); 
r = [5 * ones(1, length(t)+1);
     2 * ones(1, length(t)+1)];
u = r;
for i = 1:length(t)
    x(:, i+1) = A_cl * x(:, i) + B_cl * u(:, i);
    y(:, i) = C_cl * x(:, i) + D_cl * u(:, i); 
end
figure('Name', 'Iterative Debug', 'NumberTitle', 'off', 'Position', [500, 450, 400, 300]);
plot(t, y(1,:));  % First output
hold on;
plot(t, y(2,:));  % Second output
hold off
xlabel('Time (s)');
ylabel('System Output');
title('Controlled CL System');
legend('z_e', 'y_e');
grid on;

% 2. Iterative
A_K = K.A;
B_K = K.B;
C_K = K.C;
D_K = K.D;
t = 0:timeStep:200;
x = zeros(length(A), length(t)+1);
u = zeros(length(B(1, :)), length(t));
e = zeros(length(C(:, 1)), length(t));
y = zeros(length(C(:, 1)), length(t));
xk = zeros(length(A_K), length(t)+1);
uk = y; % property of H inf
yk = zeros(length(C_K(:, 1)), length(t)); % property of H inf  
r = [5; 2]; 
for i = 2:length(t)
    e(:, i) = r - yk(:, i-1);
    % Update system
    x(:, i+1) = A * x(:, i) + B * e(:, i);
    y(:, i) = C * x(:, i) + D * e(:, i); 
    % Update controller
    xk(:, i+1) = A_K * xk(:, i) + B_K * y(:, i);
    yk(:, i) = C_K * xk(:, i) + D_K * y(:, i);
end
figure('Name', 'Iterative 1.0', 'NumberTitle', 'off', 'Position', [100, 50, 400, 300]);
plot(t, y(1,:));  % First output
hold on;
plot(t, y(2,:));  % Second output
plot(t, r(1, :)*ones(1, length(t)), 'k--');
plot(t, r(2, :)*ones(1, length(t)), 'k--');
hold off
xlabel('Time (s)');
ylabel('System Output');
title('Controlled CL System');
legend('z_e', 'y_e');
grid on;

% 3. Iterative Alternative
A_K = K.A;
B_K = K.B;
C_K = K.C;
D_K = K.D;
t = 0:timeStep:200;
x = zeros(length(A), length(t)+1);
u = zeros(length(B(1, :)), length(t));
e = zeros(length(C(:, 1)), length(t));
y = zeros(length(C(:, 1)), length(t));
xk = zeros(length(A_K), length(t)+1);
uk = y; % property of H inf
yk = zeros(length(C_K(:, 1)), length(t)); % property of H inf  
r = [5; 2]; 
for i = 2:length(t)
    % Update controller
    xk(:, i+1) = A_K * xk(:, i) + B_K * y(:, i-1);
    yk(:, i) = C_K * xk(:, i) + D_K * y(:, i-1);
    % Get error 
    e(:, i) = r - yk(:, i);
    % Update system
    x(:, i+1) = A * x(:, i) + B * e(:, i);
    y(:, i) = C * x(:, i) + D * e(:, i); 
end
figure('Name', 'Iterative 2.0', 'NumberTitle', 'off', 'Position', [500, 50, 400, 300]);
plot(t, y(1,:));  % First output
hold on;
plot(t, y(2,:));  % Second output
plot(t, r(1, :)*ones(1, length(t)), 'k--');
plot(t, r(2, :)*ones(1, length(t)), 'k--');
hold off
xlabel('Time (s)');
ylabel('System Output');
title('Controlled CL System');
legend('z_e', 'y_e');
grid on;