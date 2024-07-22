%% This file is created to learn the different cooridnate
clear
close all

%% Basic information definition
U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;     % Rotor diameter
simTime = 6000;     % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds
Str = 0.3;                          % Strouhal number
Helix_amplitude = 6;                % Helix amplitude                
Freq = Str*U_inflow/D_IEA15MW;      % From Str, in Hz
omega_e = Freq*2*pi;

% Tranform matrix to Helix coordinate frame
% R_helix = [1 0 0; 
%            0 cos(Freq) -sin(Freq); 
%            0 sin(Freq) cos(Freq)];
% inverseR_helix = [1 0 0; 
%                   0 cos(Freq) sin(Freq); 
%                   0 -sin(Freq) cos(Freq)];

%% From Fixed Frame --> Helix Frame --> Fixed Frame, CW

% Genereate the tilt and yaw signal (fixed frame)
t = linspace(timeStep, simLen, simTime);
cutPoint = simLen / 3;
amplitudeTilt = Helix_amplitude*(t<=cutPoint) + 1.5*Helix_amplitude*(t>cutPoint);

sigTilt = Helix_amplitude * sin(omega_e*t);          
sigYaw = Helix_amplitude * sin(omega_e*t + pi/2);  % CW

% Transfer to helix frame
thetaTilt_helixFrame_store = zeros(simTime, 1);
thetaYaw_helixFrame_store = zeros(simTime, 1);
thetaTilt_fixFrame_store = zeros(simTime, 1);
thetaYaw_fixFrame_store = zeros(simTime, 1);

for i = 1:1:simTime
    R_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
               sin(omega_e*t(i)) cos(omega_e*t(i))];
    invR_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
                  -sin(omega_e*t(i)) cos(omega_e*t(i))];
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);

    % Fixed Frame ---> Helix Frame
    MtMy_helix = R_helix * [theta_tilt; theta_yaw]; 
    thetaTilt_helixFrame_store(i) = MtMy_helix(1);
    thetaYaw_helixFrame_store(i) = MtMy_helix(2);

    % Helix Frame ---> Fixed Frame
    MtMy_fix = invR_helix * MtMy_helix; 
    thetaTilt_fixFrame_store(i) = MtMy_fix(1);
    thetaYaw_fixFrame_store(i) = MtMy_fix(2);
end

figure()
subplot(3, 1, 1)
plot(t, sigTilt)
hold on
plot(t, sigYaw)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title('Fixed Frame')
legend('M_{tilt}','M_{yaw}')

subplot(3, 1, 2)
plot(t, thetaTilt_helixFrame_store)
hold on
plot(t, thetaYaw_helixFrame_store)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title("Helix Frame")
legend('M^e_{tilt}','M^e_{yaw}')

subplot(3, 1, 3)
plot(t, sigTilt)
hold on
plot(t, sigYaw)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title("Again Fixed Frame")
legend('M_{tilt}','M_{yaw}')

%% From Fixed Frame --> Helix Frame --> Fixed Frame, CCW

% Genereate the tilt and yaw signal (fixed frame)
t = linspace(timeStep, simLen, simTime);
cutPoint = simLen / 3;
amplitudeTilt = Helix_amplitude*(t<=cutPoint) + 1.5*Helix_amplitude*(t>cutPoint);

sigTilt = Helix_amplitude * sin(omega_e*t);          
sigYaw = Helix_amplitude * sin(omega_e*t - pi/2);  % CCW

% Transfer to helix frame
thetaTilt_helixFrame_store = zeros(simTime, 1);
thetaYaw_helixFrame_store = zeros(simTime, 1);
thetaTilt_fixFrame_store = zeros(simTime, 1);
thetaYaw_fixFrame_store = zeros(simTime, 1);

for i = 1:1:simTime
    R_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
               -sin(omega_e*t(i)) cos(omega_e*t(i))];
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);

    % Fixed Frame ---> Helix Frame
    MtMy_helix = R_helix * [theta_tilt; theta_yaw]; 
    thetaTilt_helixFrame_store(i) = MtMy_helix(1);
    thetaYaw_helixFrame_store(i) = MtMy_helix(2);

    % Helix Frame ---> Fixed Frame
    MtMy_fix = invR_helix * MtMy_helix; 
    thetaTilt_fixFrame_store(i) = MtMy_fix(1);
    thetaYaw_fixFrame_store(i) = MtMy_fix(2);
end

figure()
subplot(3, 1, 1)
plot(t, sigTilt)
hold on
plot(t, sigYaw)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title('Fixed Frame')
legend('M_{tilt}','M_{yaw}')

subplot(3, 1, 2)
plot(t, thetaTilt_helixFrame_store)
hold on
plot(t, thetaYaw_helixFrame_store)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title("Helix Frame")
legend('M^e_{tilt}','M^e_{yaw}')

subplot(3, 1, 3)
plot(t, sigTilt)
hold on
plot(t, sigYaw)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title("Again Fixed Frame")
legend('M_{tilt}','M_{yaw}')

%% Results of Change: Fixed Frame --> Helix Frame, CCW

% Genereate the tilt and yaw signal (fixed frame)
t = linspace(timeStep, simLen, simTime);
cutPoint1 = simLen / 3;
cutPoint2 = simLen * 2 / 3;
amplitudeChange = Helix_amplitude*(t<=cutPoint1)+2*Helix_amplitude*(t>cutPoint1);
amplitudeChange2 = Helix_amplitude*(t<=cutPoint1)+4*Helix_amplitude*(t>cutPoint1);

amplitudeChange = Helix_amplitude * (1 - t / simLen);  % Linear ramp
amplitudeChange2 = Helix_amplitude * (1 - t / simLen);  % Scaled linear ramp

% sigTilt = Helix_amplitude * sin(omega_e*t);
% sigYaw = Helix_amplitude * sin(omega_e*t - pi/2);  % CCW
sigTilt = Helix_amplitude .* sin(omega_e*t);
sigYaw = amplitudeChange .* sin(omega_e*t - pi/2);  % CCW

% Transfer to helix frame
thetaTilt_helixFrame_store = zeros(simTime, 1);
thetaYaw_helixFrame_store = zeros(simTime, 1);

for i = 1:1:simTime
    R_helix = [cos(omega_e*t(i)) sin(omega_e*t(i)); 
               -sin(omega_e*t(i)) cos(omega_e*t(i))];
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);

    % Fixed Frame ---> Helix Frame
    MtMy_helix = R_helix * [theta_tilt; theta_yaw]; 
    thetaTilt_helixFrame_store(i) = MtMy_helix(1);
    thetaYaw_helixFrame_store(i) = MtMy_helix(2);
end

figure()
subplot(2, 1, 1)
plot(t, sigTilt)
hold on
plot(t, sigYaw)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title('Fixed Frame')
legend('M_{tilt}','M_{yaw}')

subplot(2, 1, 2)
plot(t, thetaTilt_helixFrame_store)
hold on
plot(t, thetaYaw_helixFrame_store)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title("Helix Frame")
legend('M^e_{tilt}','M^e_{yaw}')

%% Results of Change: Helix Frame ---> Fixed Frame, CCW

% Genereate the tilt and yaw signal (fixed frame)
t = linspace(timeStep, simLen, simTime);
cutPoint1 = simLen / 3;
amplitudeChange = Helix_amplitude*(t<=cutPoint1)+2*Helix_amplitude*(t>cutPoint1);

sigTilt_e = 0 * ones(simTime, 1);
% sigTilt_e = [linspace(0, 4, simTime*9/20) linspace(4, 0, simTime*9/20) 0*ones(1, simTime/10)];
% sigTilt_e = [0*ones(1, simTime/5) 1*ones(1, simTime/5) 2*ones(1, simTime/5) 3*ones(1, simTime/5) 4*ones(1, simTime/5)];
% sigYaw_e = -2 * ones(simTime, 1);
sigYaw_e = [-2*ones(1, simTime/10) linspace(-2, 2, simTime*4/5) 2*ones(1, simTime/10)];

% sigTilt_e = [0*ones(1, simTime/6) 1*ones(1, simTime/6) 2*ones(1, simTime/6) 3*ones(1, simTime/6) 4*ones(1, simTime/3)];
% sigYaw_e = -2 * ones(simTime, 1);

% Transfer to helix frame
thetaTilt_fixFrame_store = zeros(simTime, 1);
thetaYaw_fixFrame_store = zeros(simTime, 1);

for i = 1:1:simTime
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    theta_tilt_e = sigTilt_e(i);
    theta_yaw_e = sigYaw_e(i);

    % Fixed Frame ---> Helix Frame
    MtMy_fixed = invR_helix * [theta_tilt_e; theta_yaw_e]; 
    thetaTilt_fixFrame_store(i) = MtMy_fixed(1);
    thetaYaw_fixFrame_store(i) = MtMy_fixed(2);
end

figure()
subplot(2, 1, 1)
plot(t, sigTilt_e)
hold on
plot(t, sigYaw_e)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title('Helix Frame')
legend('M^e_{tilt}','M^e_{yaw}')

subplot(2, 1, 2)
plot(t, thetaTilt_fixFrame_store)
hold on
plot(t, thetaYaw_fixFrame_store)
hold off
xlabel('Time [s]')
ylabel('Magnitude')
title("Fixed Frame")
legend('M_{tilt}','M_{yaw}')
