addpath('.\Functions');

%%
% Ts_prbn = 0.1;
trainData = 'trainingSet.mat';   % train set
testData = 'testingSet.mat';   % test set
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
u_train = u_train(601:end, :);
y_train = y_train(601:end, :);
u_test = u_test(601:end, :);
y_test = y_test(601:end, :);

% Time shift the signal
DeadtimeDelay = 126;

u_train = u_train(1:end-DeadtimeDelay, :);
y_train = y_train(DeadtimeDelay+1:end, :);
N_train = length(u_train);
t_train = (0:N_train-1);

u_test = u_test(1:end-DeadtimeDelay, :);
y_test = y_test(DeadtimeDelay+1:end, :);
N_test = length(u_test);
t_test = (0:N_test-1);

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

[Ga,ws] = spa_avf(u_train,y_train,1,25,[],[],'hamming');
Ga = frd(Ga,ws);
figure 
bodemag(Ga,'m');
legend('SPA AVF');
% bodemag(OLi(1:2,1:2), Ga, 'm');
% legend('real', 'spa avf');

%% System IDE in open loop 
% PBSID-varx
n = 4;
f = 30;
p = 30;

[S,X] = dordvarx(u_train,y_train,f,p,'tikh','gcv');
figure, semilogy(S,'*');
x = dmodx(X,n);
[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,u_train,y_train,f,p);

% PBSID-varmax


%% SystemIDE Results Verification
% trainset
OLi = ss(Ai,Bi,Ci,Di,1);
yi = lsim(OLi,u_train,t_train); 
disp('VAF with PBSID-varx (open loop)')
vaf(y_train, yi)                

% testset
yi_test = lsim(OLi,u_test,t_test); 
disp('VAF with PBSID-varx (open loop)')
vaf(y_test, yi_test)                