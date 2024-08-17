clear
close all
addpath('.\Functions');

%%
trainData = 'train_30min_2bw.mat';   % train set
testData = 'test_2steps.mat';        % test set
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\sysIDE\';
IDEdata_train = load([turbineName caseName trainData]);
IDEdata_test = load([turbineName caseName testData]);

%% Get Training data and Testing data
u_train = IDEdata_train.HF_beta;
y_train = IDEdata_train.HF_helixCenter_filtered;
u_test = IDEdata_test.HF_beta;
y_test = IDEdata_test.HF_helixCenter_filtered;

% Remove first few data
shiftNum = 1000;
u_train = u_train(shiftNum:end, :);
y_train = y_train(shiftNum:end, :);
u_test = u_test(shiftNum:end, :);
y_test = y_test(shiftNum:end, :);

% Time shift the signal
DeadtimeDelay = 110;
u_train = u_train(1:end-DeadtimeDelay, :);
y_train = y_train(DeadtimeDelay+1:end, :);
N_train = length(u_train);
t_train = (0:N_train-1);

u_test = u_test(1:end-DeadtimeDelay, :);
y_test = y_test(DeadtimeDelay+1:end, :);
N_test = length(u_test);
t_test = (0:N_test-1);

% Detrend data
u_train = u_train - mean(u_train);  % [data2_d,Tr] = detrend(data2);
y_train = y_train - mean(y_train);
u_test = u_test - mean(u_test);
y_test = y_test - mean(y_test);

%% 
% % Power Spectrum Density
% [M1,F1] = pwelch(u_train(:, 1),[],[],[],1/Ts_prbn);
% [M2,F2] = pwelch(u_train(:, 2),[],[],[],1/Ts_prbn);
% figure
% semilogx(F1,mag2db(M1),'k','LineWidth',1)
% hold on
% semilogx(F2,mag2db(M2),'r','LineWidth',1)
% hold off
% xlabel('Frequency [Hz]');
% ylabel('Amplitude [dB]');
% legend('\beta^e_{tilt}', '\beta^e_{yaw}')
% title('Input PSD')
% 
% [M1,F1] = pwelch(y_train(:, 1),[],[],[],1/Ts_prbn);
% [M2,F2] = pwelch(y_train(:, 2),[],[],[],1/Ts_prbn);
% figure
% semilogx(F1,mag2db(M1),'k','LineWidth',1)
% hold on
% semilogx(F2,mag2db(M2),'r','LineWidth',1)
% hold off
% xlabel('Frequency [Hz]');
% ylabel('Amplitude [dB]');
% legend('z_f', 'y_f')
% title('Output PSD')

%% PBSID-varx
n_varx = 4;     % 4
f_varx = 30;    
p_varx = 30;

[S,X] = dordvarx(u_train,y_train,f_varx,p_varx,'tikh','gcv');
figure, semilogy(S,'*');
title('Singular Value PBSID-varx')
x = dmodx(X,n_varx);
[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,u_train,y_train,f_varx,p_varx);
% State-space model
OLi = ss(Ai,Bi,Ci,Di,1);
OLi.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLi.OutputName = {'z_e','y_e'};

%% PBSID-varmax
n_varmax = 4;
f_varmax = 30;
p_varmax = 30;

[S,X] = dordvarmax(u_train,y_train,f_varmax,p_varmax,'els',1e-6,'tikh','gcv');
figure, semilogy(S,'*');
title('Singular Value PBSID-varmax')
x = dmodx(X,n_varmax);
[Av,Bv,Cv,Dv,Kv] = dx2abcdk(x,u_train,y_train,f_varmax,p_varmax);
% State-space model
OLv = ss(Av,Bv,Cv,Dv,1);
OLv.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLv.OutputName = {'z_e','y_e'};

%% SystemIDE Results Verification
% ====== trainset 
% frequency respond
[Ga,ws] = spa_avf(u_train,y_train,1,25,[],[],'hamming');
Ga = frd(Ga,ws);
figure 
bodemag(OLi, 'c', OLv, 'g', Ga, 'm');
legend('PBSID-varx', 'PBSID-varmax', 'Real spa avf');
% VAF
yi = lsim(OLi,u_train,t_train); 
yv = lsim(OLv,u_train,t_train);
disp('[Training] VAF with PBSID-varx (open loop)')
vaf(y_train, yi)            
disp('[Training] VAF with PBSID-varmax (open loop)')
vaf(y_train, yv)  

% ====== testset
% % frequency respond
% [Ga,ws] = spa_avf(u_test,y_test,1,25,[],[],'hamming');
% Ga = frd(Ga,ws);
% figure 
% bodemag(OLi, 'c', OLv, 'g', Ga, 'm');
% legend('PBSID-varx', 'PBSID-varmax', 'Real spa avf');
% VAF
yi_test = lsim(OLi,u_test,t_test); 
yv_test = lsim(OLv,u_test,t_test);
disp('================================================')
disp('[Testing] VAF with PBSID-varx (open loop)')
vaf(y_test, yi_test)         
disp('[Testing] VAF with PBSID-varmax (open loop)')
vaf(y_test, yv_test)   