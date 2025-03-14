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

%% Basic system property
fprintf('======== System property \n');
fprintf(' System Dimension: %.0f \n', size(A, 1));
fprintf(' Ctrb Matrix Rank: %.0f \n', rank(ctrb(A, B)));
fprintf(' Obsv Matrix Rank: %.0f \n', rank(obsv(A, C)));
fprintf(' Eigenvalues of A: \n');
disp(eig(A));

figure('Name', 'Bode Diagram OL System', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
bode(G);
bw = calculateBandwidth(G(1,1));   % Hz
bw_rad = bw*2*pi;
grid on
title('Bode Diagram of Open-Loop System')

%% PID Controller Design
% func: pidtune
wc_hz = 0.5*bw;
wc = wc_hz * 2*pi; % rad/timeUnit omega = 2pi*f
printOption = 'N';
bdOption = 'N';
[Ci1, info1] = pidtune(G(1,1), 'I', wc);
[Ci2, info2] = pidtune(G(1,1), 'I');
[Cpi1, info3] = pidtune(G(1,1), 'PI', wc);
[Cpi2, info4] = pidtune(G(1,1), 'PI');
C_mimo = [Cpi1, 0;
          0, 0];
C_mimo2 = [0, 0;
           0, Cpi1];
C_mimo3 = [Cpi1, 0;
           0, Cpi1];


%%%
% Open loop transfer matrix becomes
%   [C11*G11 C11*G12]
%   [C22*G21 C22*G22]
%%%
OL_ctrl = C_mimo*G;
OL_ctrl2 = C_mimo2*G;
OL_ctrl3 = C_mimo3*G;

% ==== Check Controller
% Frequency Domain
% 1. Bode Diagram of Controllers
if strcmp(bdOption, 'Y') 
    figure('Name', 'BD Controllers', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bode(Ci1, Ci2, Cpi1, Cpi2)
    axesHandles = findall(gcf, 'Type', 'axes');
    set(axesHandles(2), 'YLim', [-180 10]);
    set(axesHandles(3), 'YLim', [-100 50]);
    hold on
    % Don't forget to convert to Hz when using below to show bw !!!
    xline(axesHandles(3), bw, 'k--', 'LineWidth', 2);
    xline(axesHandles(3), wc_hz, 'k-.', 'LineWidth', 2);
    xline(axesHandles(2), bw, 'k--', 'LineWidth', 2);
    xline(axesHandles(2), wc_hz, 'k-.', 'LineWidth', 2);
    hold off
    title('Controller Bode Diagram')
    legend('I1', 'I2','PI1','PI2','\omega_b','\omega_c','Location','southeast');
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
    grid on
end

% 2. Bode Diagram of OL System G11 and G22 only
% Ci1
if strcmp(printOption, 'Y')
    fprintf('======== I1 \n');
    [Gm, Pm, ~, Wcp] = margin(G(1,1)*Ci1);
    [Gm2, Pm2, ~, Wcp2] = margin(G(1,2)*Ci1);
    [Gm3, Pm3, ~, Wcp3] = margin(G(2,1)*Ci1);
    [Gm4, Pm4, ~, Wcp4] = margin(G(2,2)*Ci1);
    fprintf(' Gain Margin: %.3f, %.3f, %.3f, %.3f [dB] \n', 20*log10(Gm), 20*log10(Gm2), 20*log10(Gm3), 20*log10(Gm4));
    fprintf(' Phase Margin: %.3f, %.3f, %.3f, %.3f [deg] \n', Pm, Pm2, Pm3, Pm4);
    fprintf(' Crossover Frequncy: %.3f, %.3f, %.3f, %.3f [rad] \n', Wcp, Wcp2, Wcp3, Wcp4);
    
    % Ci2
    fprintf('======== I2 \n');
    [Gm, Pm, ~, Wcp] = margin(G(1,1)*Ci2);
    [Gm2, Pm2, ~, Wcp2] = margin(G(1,2)*Ci2);
    [Gm3, Pm3, ~, Wcp3] = margin(G(2,1)*Ci2);
    [Gm4, Pm4, ~, Wcp4] = margin(G(2,2)*Ci2);
    fprintf(' Gain Margin: %.3f, %.3f, %.3f, %.3f [dB] \n', 20*log10(Gm), 20*log10(Gm2), 20*log10(Gm3), 20*log10(Gm4));
    fprintf(' Phase Margin: %.3f, %.3f, %.3f, %.3f [deg] \n', Pm, Pm2, Pm3, Pm4);
    fprintf(' Crossover Frequncy: %.3f, %.3f, %.3f, %.3f [rad] \n', Wcp, Wcp2, Wcp3, Wcp4);
    
    % Cpi1
    fprintf('======== PI1 \n');
    [Gm, Pm, ~, Wcp] = margin(G(1,1)*Cpi1);
    [Gm2, Pm2, ~, Wcp2] = margin(G(1,2)*Cpi1);
    [Gm3, Pm3, ~, Wcp3] = margin(G(2,1)*Cpi1);
    [Gm4, Pm4, ~, Wcp4] = margin(G(2,2)*Cpi1);
    fprintf(' Gain Margin: %.3f, %.3f, %.3f, %.3f [dB] \n', 20*log10(Gm), 20*log10(Gm2), 20*log10(Gm3), 20*log10(Gm4));
    fprintf(' Phase Margin: %.3f, %.3f, %.3f, %.3f [deg] \n', Pm, Pm2, Pm3, Pm4);
    fprintf(' Crossover Frequncy: %.3f, %.3f, %.3f, %.3f [rad] \n', Wcp, Wcp2, Wcp3, Wcp4);
    
    % Cpi2
    fprintf('======== PI2 \n');
    [Gm, Pm, ~, Wcp] = margin(G(1,1)*Cpi2);
    [Gm2, Pm2, ~, Wcp2] = margin(G(1,2)*Cpi2);
    [Gm3, Pm3, ~, Wcp3] = margin(G(2,1)*Cpi2);
    [Gm4, Pm4, ~, Wcp4] = margin(G(2,2)*Cpi2);
    fprintf(' Gain Margin: %.3f, %.3f, %.3f, %.3f [dB] \n', 20*log10(Gm), 20*log10(Gm2), 20*log10(Gm3), 20*log10(Gm4));
    fprintf(' Phase Margin: %.3f, %.3f, %.3f, %.3f [deg] \n', Pm, Pm2, Pm3, Pm4);
    fprintf(' Crossover Frequncy: %.3f, %.3f, %.3f, %.3f [rad] \n', Wcp, Wcp2, Wcp3, Wcp4);
end

% Time Domain: Step Response
closed_loop_sys = feedback(OL_ctrl, eye(2));
closed_loop_sys2 = feedback(OL_ctrl2, eye(2));
closed_loop_sys.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
closed_loop_sys.OutputName = {'z_e','y_e'};
closed_loop_sys2.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
closed_loop_sys2.OutputName = {'z_e','y_e'};
t = 0:timeStep:200;  % Time vector for simulation
figure('Name', 'After Control CL Step', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
step(closed_loop_sys, t);
hold on
step(closed_loop_sys2, t);
hold off
legend('Tilt Only','Yaw Only','Location','southeast');
h = findall(gcf, 'Type', 'axes');
set(h, 'XLim', [0 200], 'YLim', [-0.5 1.5]);  
xticks(0:50:200);
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)

% Time Domain: Step Response2
closed_loop_sys3 = feedback(OL_ctrl3, eye(2));
closed_loop_sys3.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
closed_loop_sys3.OutputName = {'z_e','y_e'};
t = 0:timeStep:249;  % Time vector for simulation
figure('Name', 'After Control CL Step', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
step(closed_loop_sys3, t);
h = findall(gcf, 'Type', 'axes');
% xticks(0:50:200);
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)

%% Explane the sequential SIMO control
% C_mimo = [Cpi1, 0;
%           0, 0];
% OL_ctrl = C_mimo*G;
% C_mimo2 = [0, 0;
%            0, Cpi1];
% OL_ctrl2 = C_mimo2*G;
% C_mimo3 = [Cpi1, 0;
%            0, Cpi1];
% OL_ctrl3 = C_mimo3*G;
% 
% figure('Name', 'BD Controllers', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% bode(OL_ctrl, OL_ctrl2, OL_ctrl3)
% axesHandles = findall(gcf, 'Type', 'axes');
% legend('OL1', 'OL2','OL3','Location','southeast');
% 
% margin(OL_ctrl3(1,1))
% margin(OL_ctrl3(1, 2))
% 
% margin(OL_ctrl3(2, 2))
% margin(OL_ctrl3(2, 1))

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
% zpk(G(1,1))
% zpk(G(1,2))
% zpk(G(2,1))
% zpk(G(2,2))

%% Step simulation (Iterative way of implementing PI controller)
% % Default way
% closed_loop_sys = feedback(OL_ctrl, eye(2));
% t = 0:timeStep:200;  % Time vector for simulation
% figure('Name', 'As a whole', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% [y, tOut] = step(closed_loop_sys, t);
% title('Controlled CL System');
% plot(tOut, y(:, 2))
% hold on
% plot(tOut, y(:, 4))
% hold off
% xlabel('Time (s)');
% ylabel('System Output');
% title('Controlled CL System');
% legend('z_e', 'y_e');
% grid on;
% 
% % iterative way
% t = 0:timeStep:200;  % Time vector for simulation
% OL_sys = buf_sys.OLi;
% x = zeros(4, length(t)+1);
% u = zeros(2, length(t));
% e = zeros(2, length(t));
% y = zeros(2, length(t));
% r = [0; 1];
% Kp = 0.557;     % 0.557; 0
% Ki = 0.0327;    % 0.0327; 0.0113           
% Kp_matrix = [0 0; 
%              0 Kp];
% Ki_matrix = [0 0; 
%              0 Ki];
% %%%!!! The problem is the way PI is implemented iteratively
% for i = 2:length(t)
%     e(:, i) = r - y(:, i-1);
%     delta_u = Ki_matrix * e(:, i) * timeStep;
%     u(:, i) = Kp_matrix * (e(:, i)-e(:, i-1)) + u(:, i-1) + delta_u;
%     x(:, i+1) = OL_sys.A * x(:, i) + OL_sys.B * u(:, i);  
%     y(:, i) = OL_sys.C * x(:, i) + OL_sys.D * u(:, i);
% end
% figure('Name', 'Iterative', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
% plot(t, y(1,:));  % First output
% hold on;
% plot(t, y(2,:));  % Second output
% hold off
% xlabel('Time (s)');
% ylabel('System Output');
% title('Controlled CL System');
% legend('z_e', 'y_e');
% grid on;