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
G = tf(buf_sys.OLi);        % transfer matrix 
timeStep = 0.1;

% Basic system property
fprintf('======== System property \n');
fprintf(' System Dimension: %.0f \n', size(A, 1));
fprintf(' Ctrb Matrix Rank: %.0f \n', rank(ctrb(A, B)));
fprintf(' Obsv Matrix Rank: %.0f \n', rank(obsv(A, C)));
fprintf(' Eigenvalues of A: \n');
disp(eig(A));

%% H infinity Control Design
W_s = tf([1, 1.6], [100, 1]);  % Emphasizes performance and disturbance rejection
W_t = tf([0.01, 1], [6, 1]);   % Emphasizes robustness and noise rejection
W_u = tf([1, 1], [50, 1]);     % Emphasizes input magnitude  
% Working weight function
% W_s = tf([1, 1.6], [100, 1]);
% W_t = tf([1, 1], [1, 1]);
W_s_d = c2d(W_s, timeStep, 'tustin');
W_t_d = c2d(W_t, timeStep, 'tustin');
W_u_d = c2d(W_u, timeStep, 'tustin');
P = augw(sys, W_s_d, W_u_d, W_t_d);  % augw creates the weighted augmented plant
ncont = 2; 
nmeas = 2; 
[K,CL,gamma] = hinfsyn(P,nmeas,ncont);
sys_cl = feedback(sys, K);

showFreqOption = 'N';
showTimeOption = 'Y';

% Frequency Response
bwG11 = calculateBandwidth(G(1, 1));
bwG22 = calculateBandwidth(G(2, 2));
bw = (bwG11 + bwG22)/2; 
L_mimo = G * K;                     % open-loop tf
S_mimo = feedback(eye(2), L_mimo);  % disturbance --- output
T_mimo = feedback(L_mimo, eye(2));  % reference --- output
U_mimo = feedback(K, G);            % disturbance --- negative control input
L_mimo.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
L_mimo.OutputName = {'z_e','y_e'};

if strcmp(showFreqOption, 'Y')
    figure('Name', 'Frequence Result', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bode(L_mimo, S_mimo, T_mimo, U_mimo)
    hold on
    axesHandles = findall(gcf, 'Type', 'axes');
    % remember to switch the unit to Hz 
    xline(axesHandles(2), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(4), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(6), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(8), bw, 'k--', 'LineWidth', 1);
    
    xline(axesHandles(3), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(5), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(7), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(9), bw, 'k--', 'LineWidth', 1);
    hold off
    grid on
    title('Frequency Response of MIMO System')
    legend('L', 'S', 'T', 'U','\omega_b', 'Location','southeast')
    % setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
end

% Study the direction 
[S, wout] = sigma(G); % singular value, rad/s
condition_number = max(S) ./ min(S);
bw_rad = bw * 2*pi;
[~, idx] = min(abs(wout - bw_rad));
S(:, idx);
condition_number(idx);

% Time Response
if strcmp(showTimeOption, 'Y')
    t = 0:0.1:200;
    figure('Name', 'Step Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    step(sys, t)
    hold on
    step(sys_cl, t)
    axesHandles = findall(gcf, 'Type', 'axes');
    yline(axesHandles(5), 1, 'k--', 'LineWidth', 2);
    yline(axesHandles(4), 0, 'k--', 'LineWidth', 2);
    yline(axesHandles(3), 0, 'k--', 'LineWidth', 2);
    yline(axesHandles(2), 1, 'k--', 'LineWidth', 2);
    hold off
    legend('OL', 'CL')
    grid on
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
end

%% H infinity controller design 2 (Matrix Form)
% s = tf('s');
% wB1 = 0.0025*2*pi;
% A = 1/100;
% M = 100;
% W_s = (s/M+wB1)/(s+wB1*A);
W_s = tf([1, 1.6], [100, 1]);  % Emphasizes performance and disturbance rejection
W_t = tf([0.01, 1], [6, 1]);   % Emphasizes robustness and noise rejection
W_u = tf([1, 1], [50, 1]);     % Emphasizes input magnitude  
% W_u = eye(2);
% Working weight function
% W_s = tf([1, 1.6], [100, 1]);
% W_t = tf([1, 1], [1, 1]);
W_s_d = c2d(W_s, timeStep, 'tustin');
W_t_d = c2d(W_t, timeStep, 'tustin');
W_u_d = c2d(W_u, timeStep, 'tustin');
W_sMatrix = blkdiag(W_s_d, W_s_d);
W_tMatrix = blkdiag(W_t_d, W_t_d);
W_uMatrix = blkdiag(W_u_d, W_u_d);
P = augw(sys, W_sMatrix, W_tMatrix, W_uMatrix);  % augw creates the weighted augmented plant
ncont = 2; 
nmeas = 2; 
[K,CL,gamma] = hinfsyn(P,nmeas,ncont);
sys_cl = feedback(sys, K);

showFreqOption = 'N';
showTimeOption = 'Y';

% Frequency Response
bwG11 = calculateBandwidth(G(1, 1));
bwG22 = calculateBandwidth(G(2, 2));
bw = (bwG11 + bwG22)/2; 
L_mimo = G * K;                     % open-loop tf
S_mimo = feedback(eye(2), L_mimo);  % disturbance --- output
T_mimo = feedback(L_mimo, eye(2));  % reference --- output
U_mimo = feedback(K, G);            % disturbance --- negative control input
L_mimo.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
L_mimo.OutputName = {'z_e','y_e'};

if strcmp(showFreqOption, 'Y')
    figure('Name', 'Frequence Result', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bode(L_mimo, S_mimo, T_mimo, U_mimo)
    hold on
    axesHandles = findall(gcf, 'Type', 'axes');
    % remember to switch the unit to Hz 
    xline(axesHandles(2), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(4), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(6), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(8), bw, 'k--', 'LineWidth', 1);
    
    xline(axesHandles(3), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(5), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(7), bw, 'k--', 'LineWidth', 1);
    xline(axesHandles(9), bw, 'k--', 'LineWidth', 1);
    hold off
    grid on
    title('Frequency Response of MIMO System')
    legend('L', 'S', 'T', 'U','\omega_b', 'Location','southeast')
    % setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
end

% Study the direction 
[S, wout] = sigma(G); % singular value, rad/s
condition_number = max(S) ./ min(S);
bw_rad = bw * 2*pi;
[~, idx] = min(abs(wout - bw_rad));
S(:, idx);
condition_number(idx);

% Time Response
if strcmp(showTimeOption, 'Y')
    t = 0:0.1:200;
    figure('Name', 'Step Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    step(sys, t)
    hold on
    step(sys_cl, t)
    axesHandles = findall(gcf, 'Type', 'axes');
    yline(axesHandles(5), 1, 'k--', 'LineWidth', 2);
    yline(axesHandles(4), 0, 'k--', 'LineWidth', 2);
    yline(axesHandles(3), 0, 'k--', 'LineWidth', 2);
    yline(axesHandles(2), 1, 'k--', 'LineWidth', 2);
    hold off
    legend('OL', 'CL')
    grid on
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
end

%% Step Response Implementation: Default vs. Iterative
% close all
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