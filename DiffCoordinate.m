%% This file is created to learn the different cooridnate
clear
close all

%% Basic information definition
U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;     % Rotor diameter
simTime = 10000;     % in timestep, actual time is simTime*timestep(Q-blade define)
timeStep = 0.1;    % same with the Q-blade setting
simLen = simTime * timeStep; % seconds
Str = 0.3;                          % Strouhal number
Helix_amplitude = 4;                % Helix amplitude                
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
    R_helix = [1 0 0; 
           0 cos(omega_e*t(i)) -sin(omega_e*t(i)); 
           0 sin(omega_e*t(i)) cos(omega_e*t(i))];
    invR_helix = [1 0 0; 
                  0 cos(omega_e*t(i)) sin(omega_e*t(i)); 
                  0 -sin(omega_e*t(i)) cos(omega_e*t(i))];
    theta_tilt = sigTilt(i);
    theta_yaw = sigYaw(i);
    theta_col = 0;

    % Fixed Frame ---> Helix Frame
    MtMy_helix = R_helix * [theta_col; theta_tilt; theta_yaw]; 
    thetaTilt_helixFrame_store(i) = MtMy_helix(2);
    thetaYaw_helixFrame_store(i) = MtMy_helix(3);

    % Helix Frame ---> Fixed Frame
    MtMy_fix = invR_helix * MtMy_helix; 
    thetaTilt_fixFrame_store(i) = MtMy_fix(2);
    thetaYaw_fixFrame_store(i) = MtMy_fix(3);
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