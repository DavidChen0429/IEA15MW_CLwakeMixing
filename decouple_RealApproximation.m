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
% G = G([2, 1], :);           % Switch row because of control relation

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

% % 1. steady-state decoupling
ss_compensator = eye(2)/(G_ss) * 3; % match gain *3
OLi = series(buf_sys.OLi, ss_compensator);
OLi.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLi.OutputName = {'z_e','y_e'};
G_ss_dcpl = dcgain(OLi);
RGA_ssdcpl = G_ss_dcpl .* (inv(G_ss_dcpl))';
disp(RGA_ssdcpl)

% 2. Specific frequency decouple / pre-compensator
bwG11 = calculateBandwidth(G(1, 1));
bwG22 = calculateBandwidth(G(2, 2));
bw = (bwG11 + bwG22)/2;  
[G_mag, G_phase] = bode(buf_sys.OLi, bw);  
G_bw = freqresp(buf_sys.OLi, bw);  % Get the frequency response matrix at bw
bw_compensator = inv(G_bw);  % Decoupling matrix at bw
G_bw_real = real(bw_compensator) * 3;  % Real approximation
OLi = series(buf_sys.OLi, G_bw_real);
OLi.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLi.OutputName = {'z_e','y_e'};
G_bw_dcpl = dcgain(OLi);
RGA_bwdcpl = G_bw_dcpl .* (inv(G_bw_dcpl))';
disp(RGA_bwdcpl)

% % Real-approximation
[V_left, D_left] = eig(buf_sys.OLi.A'); 
V1 = V_left(:, 1);
V2 = V_left(:, 2);
V3 = V_left(:, 3);
V4 = V_left(:, 4);
V = [V1';V2';V3';V4'];
Alpha1 = real(V1);
Beta1 = imag(V1);
Alpha2 = real(V2);
Beta2 = imag(V2);
Alpha3 = real(V3);
Beta3 = imag(V3);
Alpha4 = real(V4);
Beta4 = imag(V4);
D1 = (Alpha2*Alpha2' + Beta2*Beta2') + (Alpha3*Alpha3' + Beta3*Beta3') + (Alpha4*Alpha4' + Beta4*Beta4');
D2 = (Alpha1*Alpha1' + Beta1*Beta1') + (Alpha3*Alpha3' + Beta3*Beta3') + (Alpha4*Alpha4' + Beta4*Beta4');
D3 = (Alpha2*Alpha2' + Beta2*Beta2') + (Alpha1*Alpha1' + Beta1*Beta1') + (Alpha4*Alpha4' + Beta4*Beta4');
D4 = (Alpha2*Alpha2' + Beta2*Beta2') + (Alpha3*Alpha3' + Beta3*Beta3') + (Alpha1*Alpha1' + Beta1*Beta1');
[R1, J1] = eig(D1);
[R2, J2] = eig(D2);
[R3, J3] = eig(D3);
[R4, J4] = eig(D4);

% % 3. dynamic decoupling
% % This currently has the issue of not being a non-minimal phase system
% % which introduces unstable RHP-poles when inverse G
% tzero(G)
% dyn_compensator = inv(G);
% decouple_sys2 = series(dyn_compensator, buf_sys.OLi);
% G_ss2 = dcgain(decouple_sys2);
% RGA2 = G_ss2 .* (inv(G_ss2))';

%% Compare result
yi2_train = lsim(buf_sys.OLi,us,t_train); % testing set
yi2d_train = lsim(OLi,us,t_train); % ss decouple 
yi2_test = lsim(buf_sys.OLi,us2,t_test); % testing set
yi2d_test = lsim(OLi,us2,t_test); % ss decouple 

% VAF
disp('=================================================')
disp('[Training] VAF with PBSID-varx (open loop)')
vaf(ys, yi2_train)
vaf(ys, yi2d_train)
disp('[Testing] VAF with PBSID-varx (open loop)')
vaf(ys2, yi2_test)  
vaf(ys2, yi2d_test)  

% Validation (Frequency domain)
[Ga,ws] = spa_avf(us,ys,timeStep,6,[],[],'hamming');
Ga = frd(Ga,ws);
figure('Name', 'Frequence Response', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
bodemag(buf_sys.OLi, OLi, Ga);
hold on
axesHandles = findall(gcf, 'Type', 'axes');
for k = 1:length(axesHandles)
    yline(axesHandles(k), 0, 'k--', 'LineWidth', 0.5); % '--r' makes it a dashed red line
end
hold off;
grid on
legend('Original Sys','Decoupled Sys', 'Real');

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
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
legend('z_e - tilt','y_e - yaw','z_{e,d} - tilt','y_{e,d} - yaw')
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