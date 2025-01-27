%% Explore an easier option of controller design
clear
close all
controllerOption = 'LL'; % LL, PID

% Infinite gain margin 
G = tf(1, [1 10 1]);

%% Alternative easier controller
% P, I, D
Tf = 0.1;
Kp = 1;
Ki = 1;
Kd = 1;
P = pid(Kp,0,0,Tf);
I = pid(0,Ki,0,Tf);
D = pid(0,0,Kd,Tf);
PI = pid(Kp,Ki,0,Tf);
PD = pid(Kp,0,Kd,Tf);
PID = pid(Kp,Ki,Kd,Tf);
% Lead-Lag 
alpha = 0.5;
Tlead = 0.1;
beta = 2;
Tlag = 0.1;
Lead = tf([Tlead, 1], [alpha*Tlead, 1]);
Lag = tf([beta*Tlag, beta], [beta*Tlag, 1]);
LeadLag = Lead*Lag;
if strcmp(controllerOption, 'PI')
    figure('Name', 'PID Frequency', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
    bode(P)
    hold on
    bode(I)
    bode(D)
    bode(PI)
    bode(PD)
    bode(PID)
    hold off
    legend('P','I','D','PI','PD','PID')
    title('PID')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
elseif strcmp(controllerOption, 'LL')
    figure('Name', 'Lead-Lag Frequency', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600])
    bode(Lead)
    hold on
    bode(Lag)
    bode(LeadLag)
    bode(Lead * PI)
    hold off
    legend('Lead','Lag','Lead-Lag','Custom')
    title('Lead-Lag')
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',2)
end