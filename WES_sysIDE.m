clear
close all
addpath('.\Functions');

%% Load model
buf_sys = load('Model/RightTransform_Azimuth96/ModelOrder4_noise1p_opposite.mat');
% decouple_sys = load('Model\ModelOrder4_AzimuthOffset.mat');
bw = 0.0145 * 2 *pi;    % rad/s
lw = 2; % line width

A = buf_sys.OLi.A;
B = buf_sys.OLi.B;
C = buf_sys.OLi.C;
D = buf_sys.OLi.D;

sys = buf_sys.OLi;
G = tf(buf_sys.OLi);        % transfer matrix 

%% Signals for system IDE
simTime = 72000*2;   % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds

% ======== Train data
bw_tilt = 0.0175;       % estimated bandwidth [Hz]
bw_yaw = 0.0175;        % estimated bandwidth [Hz]
N_signal = simTime;     % signal length [s] simLen
AMPL_signal = 1;        % amplitude
% === Pseudoradom Binary
% Beta_tilt
N_prbn = simTime;       % signal length [s] simLen
AMPL_prbn = 1;          % amplitude
Ts_prbn = timeStep;     % sampling time [s] timeStep
F_prbn = 1*bw_tilt;     % cutoff frequency [Hz]
Fstop_prbn = 1*bw_tilt; % band-stop filtered around this frequency
T0_prbn = 0;            % starting time [s]
P_prbn = 2;             % number of channels
IDEsig_tilt = idprbs(N_prbn,AMPL_prbn,Ts_prbn,F_prbn,Fstop_prbn,T0_prbn,P_prbn);
ns_prbn = floor((length(IDEsig_tilt)-N_prbn)/2);
sigTilt_e = IDEsig_tilt(ns_prbn+1:N_prbn+ns_prbn,1);   % tailor length
sigYaw_e = IDEsig_tilt(ns_prbn+1:N_prbn+ns_prbn,2);    % tailor length
% Beta_yaw
% F_prbn = 1*bw_yaw;      % cutoff frequency [Hz]
% Fstop_prbn = 1*bw_yaw;  % band-stop filtered around this frequency
% T0_prbn = 0;            % starting time [s]
% P_prbn = 1;             % number of channels
% IDEsig_yaw = idprbs(N_prbn,AMPL_prbn,Ts_prbn,F_prbn,Fstop_prbn,T0_prbn,P_prbn);
% ns_prbn = floor((length(IDEsig_yaw)-N_prbn)/2);
% sigYaw_e = IDEsig_yaw(ns_prbn+1:N_prbn+ns_prbn,1);    % tailor length

% === Chirp
% signal_length = simTime;      
% t = 0:timeStep:simLen;  
% f0 = 0.5*bw;                   
% f1 = 1.5*bw;                  
% sigTilt_e = chirp(t, f0, simLen, f1) * AMPL_signal;
% sigYaw_e = chirp(t, f0, simLen, f1) * AMPL_signal;

% === Add disturbances (Gaussian noise)
disturbance = randn(N_signal, 2);                   % noise
noise_level = AMPL_signal * 0.03; % Adjust the noise level as needed
noise_tilt = noise_level * randn(size(sigTilt_e));
noise_yaw = noise_level * randn(size(sigYaw_e));
sigTilt_e = sigTilt_e + noise_tilt;
sigYaw_e = sigYaw_e + noise_yaw;

% [u1s,Du1,u2s,Du2] = sigscale(sigTilt_e,sigYaw_e); % signal scaling

% % ======== Test data
% Helix_amplitude = 1;
% steps = [0*ones(1, simTime/5) Helix_amplitude*ones(1, simTime/5) 0*ones(1, simTime/5) Helix_amplitude*ones(1, simTime/5) 0*ones(1, simTime/5)];
% sigTilt_e = steps;                  % 0 * ones(simTime, 1)
% sigYaw_e = steps;                   % 0 * ones(simTime, 1)
% step_size = 1000;
% num_steps = simTime / step_size;
% Helix_amplitude = 1;
% steps = repmat([0 Helix_amplitude], 1, num_steps/2);
% original_time = linspace(0, 1, length(steps));
% new_time = linspace(0, 1, simTime);
% stretched_signal = interp1(original_time, steps, new_time, 'previous');
% steps = stretched_signal;
% sigTilt_e = steps;   % 0 * ones(simTime, 1)
% sigYaw_e = 0 * ones(simTime, 1);    % 0 * ones(simTime, 1)

% Power Spectrum Density
option = 'N';
if strcmp(option, 'Y')
    [M1,F1] = pwelch(sigTilt_e,[],[],[],1/timeStep);
    [M2,F2] = pwelch(sigYaw_e,[],[],[],1/timeStep);
    figure
    semilogx(F1,mag2db(M1),'k','LineWidth',1)
    hold on
    semilogx(F2,mag2db(M2),'r','LineWidth',1)
    yline(0, '--', 'LineWidth', 1)
    hold off
    xlabel('Frequency [Hz]');
    ylabel('Amplitude [dB]');
    legend('\beta^e_{tilt}', '\beta^e_{yaw}')
    title('Input PSD')
end

% Calculate the NSR (Noise-Signal-Ratio)
noise_tilt_power = mean(noise_tilt.^2);
noise_yaw_power = mean(noise_yaw.^2);
signal_tilt_power = mean((sigTilt_e - noise_tilt).^2);
signal_yaw_power = mean((sigYaw_e - noise_yaw).^2);
NSR_tilt = noise_tilt_power / signal_tilt_power;
NSR_yaw = noise_yaw_power / signal_yaw_power;
disp(['Noise-to-Signal Ratio Tilt: ', num2str(NSR_tilt)]);
disp(['Noise-to-Signal Ratio Yaw: ', num2str(NSR_yaw)]);

%% Get the output of G
us = [sigTilt_e, sigYaw_e];  % concatenate signals
% Simulate the system response
T = (0:simTime-1) * timeStep;
[ys, t] = lsim(G, us, T);
noise_level2 = 0.03;
ys(:, 1) = ys(:, 1) + noise_level2 * randn(size(ys(:, 1))); % z
ys(:, 2) = ys(:, 2) + noise_level2 * randn(size(ys(:, 2))); % y

%% Bode Diagram
[Ga,ws] = spa_avf(us,ys,timeStep,8,[],[],'hamming');
Ga = frd(Ga,ws);
w = logspace(-3, 0, 500);  % from 0.01 to 1 rad/s

figure('Name', 'Frequence Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
bodemag(Ga, G, w)
hold on
axesHandles = findall(gcf, 'Type', 'axes');
% yline(axesHandles(3), 0, 'k', 'LineWidth', lw);
% yline(axesHandles(5), 0, 'k', 'LineWidth', lw);
% yline(axesHandles(7), 0, 'k', 'LineWidth', lw);
% yline(axesHandles(9), 0, 'k', 'LineWidth', lw);
xline(axesHandles(3), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(5), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(7), bw, 'k--', 'LineWidth', lw);
xline(axesHandles(9), bw, 'k--', 'LineWidth', lw);
hold off;
grid on;
legend('Real','Identified Sys','f_{bw}','Location','southwest');
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',lw)