% Demo to show non-parametric identification methods in closed-loop
clear, clc, close all;

% parameters
N = 4000; % number of samples
n = 4;    % order of system

% state-space matrices
A = [0.67 0.67 0 0; -0.67 0.67 0 0; 0 0 -0.67 -0.67; 0 0 0.67 -0.67];
B = [0.6598 -0.5256 -0.6968 -0.1474
    1.9698 0.4845 0.1722 0.5646
    4.3171 -0.4879 0.6484 -0.4660
    -2.6436 -0.3416 -0.9400 0.1032]';
C = [-0.3749 0.0751 -0.5225 0.5830; -0.8977 0.7543 0.1159 0.0982];
D = [zeros(2) eye(2)];

% open-loop system
OL = ss(A,B,C,D,1);

% closed-loop system
F = diag([0.25 0.25]);
CL = feedback(OL,F,[1 2],[1 2],-1);

% input signals
t = 0:N-1;      % time samples
r = randn(2,N); % reference signal
u = randn(2,N); % input signal

% noise
SNR = 20;  % the signal to noise ratio (var(y) = 1)
vare = sqrt(10^(-SNR/20));
e = vare.*randn(2,N); % noise signal


% simulation of open loop
x0 = zeros(n,1);
y = lsim(OL,[u; e]',t',x0)';

% ETFE (matlab)
dat = iddata(y',u',1);
Ge = etfe(dat,50,2048);
Ge.Notes = 'bug';
Ge = frd(Ge);

% SPA (matlab)
Gs = spa(dat,50);
Gs.Notes = 'bug';
Gs = frd(Gs);

% SPA with freq. averaging (open-loop)
[Ga,ws] = spa_avf(u,y,1,25,[],[],'hamming');
Ga = frd(Ga,ws);

% simulation (open loop)
figure, bodemag(OL(1:2,1:2),Ge,'c',Gs,'g',Ga,'m');
legend('REAL','ETFE (ident)','SPA (ident)','SPA AVF');


% simulation of closed loop
y = lsim(CL,[r; e]',t',x0)';
e = (r - F*y);

% ETFE (matlab)
dat = iddata(y',e',1);
Ge = etfe(dat,50,2048);
Ge.Notes = 'Bug';
Ge = frd(Ge);

% SPA (matlab)
Gs = spa(dat,50);
Gs.Notes = 'Bug';
Gs = frd(Gs);

% SPA with freq. averaging (closed-loop)
[Ga,ws] = spa_avf(e,y,r,1,25,[],[],'hamming');
Ga = frd(Ga,ws);

% simulation (closed loop)
figure, bodemag(OL(1:2,1:2),Ge,'c',Gs,'g',Ga,'m');
legend('REAL','ETFE (ident)','SPA (ident)','SPA AVF');