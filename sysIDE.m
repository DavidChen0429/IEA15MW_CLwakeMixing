clear
close all
addpath('.\Functions');

%% Get Training data and Testing data
trainData = 'train_30min_1bw.mat';       % train set
testData = 'stepResponse.mat';            % test set
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\sysIDE\';
IDEdata_train = load([turbineName caseName trainData]);
IDEdata_test = load([turbineName caseName testData]);

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
t_train = (0:N_train-1);

u_test = u_test(1:end-DeadtimeDelay, :);
y_test = y_test(DeadtimeDelay+1:end, :);
N_test = length(u_test);
t_test = (0:N_test-1);

% Detrend data
u_train = detrend(u_train, 'constant');
y_train = detrend(y_train, 'constant');
u_test = detrend(u_test, 'constant');
y_test = detrend(y_test, 'constant');

% Signal scaling
[us,Du,ys,Dy] = sigscale(u_train, y_train);
[us2,Du2,ys2,Dy2] = sigscale(u_test, y_test);

%% Power Spectrum Density
Ts_prbn = 0.1;
[M1,F1] = pwelch(us(1, :),[],[],[],1/Ts_prbn);
[M2,F2] = pwelch(us(2, :),[],[],[],1/Ts_prbn);
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

[M1,F1] = pwelch(ys(1, :),[],[],[],1/Ts_prbn);
[M2,F2] = pwelch(ys(2, :),[],[],[],1/Ts_prbn);
figure
semilogx(F1,mag2db(M1),'k','LineWidth',1)
hold on
semilogx(F2,mag2db(M2),'r','LineWidth',1)
yline(0, '--', 'LineWidth', 1)
hold off
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');
legend('z_f', 'y_f')
title('Output PSD')

%% PBSID-varx
n_varx = 10;
f_varx = 30;    
p_varx = 30;

[S,X] = dordvarx(us,ys,f_varx,p_varx,'tikh','gcv');
% figure, semilogy(S,'*');
% title('Singular Value PBSID-varx')
x = dmodx(X,n_varx);

% system IDE
[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,us,ys,f_varx,p_varx);
OLi = ss(Ai,Bi,Ci,Di,1);
OLi.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLi.OutputName = {'z_e','y_e'};

% Validation (Frequency domain)
[Ga,ws] = spa_avf(us,ys,1,25,[],[],'hamming');
Ga = frd(Ga,ws);
figure 
bodemag(OLi, 'c', Ga, 'm');
legend('PBSID-varx', 'SPA');
yi = lsim(OLi,us,t_train); 
disp('[Training] VAF with PBSID-varx (open loop)')
vaf(ys, yi)   

% yi_test = lsim(OLi,u_test,t_test);
% disp('[Testing] VAF with PBSID-varx (open loop)')
% vaf(y_test, yi_test)  

% Validation (Time domain)
yi2 = lsim(OLi,us,t_train); % training set
figure()
subplot(2,1,1)
plot(yi2(:, 1))
hold on
plot(ys(1, :))
hold off
legend('predict', 'real')
subplot(2,1,2)
plot(yi2(:, 2))
hold on
plot(ys(2, :))
hold off
legend('predict', 'real')

yi2 = lsim(OLi,us2,t_test); % testing set
figure()
subplot(2,1,1)
plot(yi2(:, 1))
hold on
plot(ys2(1, :))
hold off
legend('predict', 'real')
subplot(2,1,2)
plot(yi2(:, 2))
hold on
plot(ys2(2, :))
hold off
legend('predict', 'real')

%% % PBSID-opt
n_opt = 10;
f_opt = 30;
p_opt = 30;

[S,x] = dordvarx(us,ys,f_opt,p_opt,'tikh','gcv');
% figure, semilogy(S,'x');
% title('Singular values')
x = dmodx(x,n_opt);

[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,us,ys,f_opt,p_opt);
%[Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,us,ys,f,p,'stable'); % forces stability
Dat = iddata(ys',us',n_opt);
Mi = abcdk2idss(Dat,Ai,Bi,Ci,Di,Ki);
set(Mi,'SSParameterization','Free','DisturbanceModel','Estimate','nk',zeros(1,2));
Mp = pem(Dat,Mi);
OLi = ss(Mi);   % prediction error method
% OLi2 = Dy*OLi*inv(Du);   % Y=DY*YS U=DU*US
OLp = ss(Mp);   % prediction error method optimization
% OLp2 = Dy*OLp*inv(Du);   % Y=DY*YS U=DU*US
OLi.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLi.OutputName = {'z_e','y_e'};
OLp.InputName = {'\beta^e_{tilt}', '\beta^e_{yaw}'};
OLp.OutputName = {'z_e','y_e'};

% Variance-accounted-for (by Kalman filter)
yest = predict(Mi,Dat);
x0 = findstates(Mi,Dat);
disp('VAF of identified system')
vaf(ys,yest.y)

yest = predict(Mp,Dat);
x0 = findstates(Mp,Dat);
disp('VAF of optimized system')
vaf(ys,yest.y)

% Frequency response
[Ga,ws] = spa_avf(us,ys,1,25,[],[],'hamming');
OLa = frd(Ga,ws);
figure, bodemag(OLi,'c', OLp, 'g', OLa,'m');
legend('PBSID-opt', 'PEM', 'SPA');   

% Validation (Time domain)
yi2 = lsim(OLi,us,t_train);
figure()
subplot(2,1,1)
plot(yi2(:, 1))
hold on
plot(ys(1, :))
hold off
legend('predict', 'real')
subplot(2,1,2)
plot(yi2(:, 2))
hold on
plot(ys(2, :))
hold off
legend('predict', 'real')

yi2 = lsim(OLi,us2,t_test);
figure()
subplot(2,1,1)
plot(yi2(:, 1))
hold on
plot(ys2(1, :))
hold off
legend('predict', 'real')
subplot(2,1,2)
plot(yi2(:, 2))
hold on
plot(ys2(2, :))
hold off
legend('predict', 'real')