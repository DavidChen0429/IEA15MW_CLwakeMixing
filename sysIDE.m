addpath('.\Functions');

%%
Fs = 10;  % sampling frequency Hz
Fc = 0.05;  % cutoff frequency Hz
fileName = 'HF_Uni_sysIDE.mat';   % Fixed Frame
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\sysIDE\';
IDEdata = load([turbineName caseName fileName]);

%% 
% Power Spectrum Density
[M1,F1] = pwelch(IDEdata.HF_beta(:, 1),[],[],[],1/Ts_prbn);
[M2,F2] = pwelch(IDEdata.HF_beta(:, 2),[],[],[],1/Ts_prbn);
figure
semilogx(F1,mag2db(M1),'k','LineWidth',1)
hold on
semilogx(F2,mag2db(M2),'r','LineWidth',1)
hold off
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');
legend('\beta^e_{tilt}', '\beta^e_{yaw}')
title('Input PSD')

[M1,F1] = pwelch(IDEdata.HF_helixCenter_filtered(:, 1),[],[],[],1/Ts_prbn);
[M2,F2] = pwelch(IDEdata.HF_helixCenter_filtered(:, 2),[],[],[],1/Ts_prbn);
figure
semilogx(F1,mag2db(M1),'k','LineWidth',1)
hold on
semilogx(F2,mag2db(M2),'r','LineWidth',1)
hold off
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');
legend('z_f', 'y_f')
title('Output PSD')

%% System IDE in open loop 
u = IDEdata.HF_beta;
y = IDEdata.HF_helixCenter_filtered;

% PBSID-varx
n = 4;
f = 10;
p = 10;
[S,X] = dordvarx(u,y,f,p,'tikh','gcv');
figure, semilogy(S,'*');
x = dmodx(X,n);
[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,u,y,f,p);

% Verification results
OLi = ss(Ai,Bi,Ci,Di,1);