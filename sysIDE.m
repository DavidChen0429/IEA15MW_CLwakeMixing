addpath('.\Functions');
addpath('.\CLSPA\CLSPA');

%%
Fs = 10;  % sampling frequency Hz
Fc = 0.05;  % cutoff frequency Hz
fileName = 'HF_Uni_sysIDE.mat';   % Fixed Frame
turbineName = '.\Data\NREL5MW\';
caseName = 'Str0.3_U10_1Dd_10Hz_CCW\sysIDE\';
IDEdata = load([turbineName caseName fileName]);
% BetaCenter_Visualization(IDEdata, Fs, Fc)

HelixCenterZ = IDEdata.HF_helixCenter(:, 1);
HelixCenterY = IDEdata.HF_helixCenter(:, 2);
HF_thetaTilt = IDEdata.HF_theta(:, 1);
HF_thetaYaw = IDEdata.HF_theta(:, 2);
u = [reshape(HF_thetaTilt, 1, length(HF_thetaYaw));
    reshape(HF_thetaYaw, 1, length(HF_thetaYaw))];
y = [reshape(HelixCenterZ, 1, length(HelixCenterZ));
    reshape(HelixCenterY, 1, length(HelixCenterY))];

% Data-processing
% Time shift 
timeDelay = 126;
u_shift = u(:, 1:end-timeDelay);
y_shift = y(:, timeDelay+1:end);
% Lowpass Filter
% y(1, :) = lowpassFilter(y(1, :), Fs, Fc);
% y(2, :) = lowpassFilter(y(2, :), Fs, Fc);

% SPA with freq. averaging (open-loop)
[Ga,ws] = spa_avf(u,y,1,25,[],[],'hamming');
Ga = frd(Ga,ws);

% simulation (open loop)
figure, bodemag(Ga,'m');
legend('SPA AVF');

%% 
data = iddata(y, u, Fs); 
na = 2; % Order of the autoregressive part
nb = 2; % Order of the exogenous input part
nk = 1; % Input-output delay
model = arx(data, [na nb nk]);


% % simulation of closed loop
% y = lsim(CL,[r; e]',t',x0)';
% e = (r - F*y);
% 
% % ETFE (matlab)
% dat = iddata(y',e',1);
% Ge = etfe(dat,50,2048);
% Ge.Notes = 'Bug';
% Ge = frd(Ge);
% 
% % SPA (matlab)
% Gs = spa(dat,50);
% Gs.Notes = 'Bug';
% Gs = frd(Gs);
% 
% % SPA with freq. averaging (closed-loop)
% [Ga,ws] = spa_avf(e,y,r,1,25,[],[],'hamming');
% Ga = frd(Ga,ws);
% 
% % simulation (closed loop)
% figure, bodemag(OL(1:2,1:2),Ge,'c',Gs,'g',Ga,'m');
% legend('REAL','ETFE (ident)','SPA (ident)','SPA AVF');