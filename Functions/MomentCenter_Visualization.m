function [] = MomentCenter_Visualization(SimData, Fs, Fc)

% Basic helix information
U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;
Str = 0.3;
Freq = Str*U_inflow/D_IEA15MW;      % From Str, in Hz
omega_e = Freq*2*pi;

LiDAR_sample= SimData.LiDAR_data;
TiltYaw_FF = SimData.thetaTiltYaw_fixedFrame;
TiltYaw_HF = SimData.thetaTiltYaw_helixFrame;

datalength = size(LiDAR_sample);
simTime = datalength(1); 
timeStep = 0.1;   
simLen = simTime * timeStep; % seconds
t = linspace(1, simTime, simTime);
t2 = linspace(1, simLen, simTime);

wakeCenterYo = arrayfun(@(x) x.centerY, LiDAR_sample);
wakeCenterZo = arrayfun(@(x) x.centerZ, LiDAR_sample);

wakeCenterYf = lowpassFilter(wakeCenterYo, Fs, Fc);
wakeCenterZf = lowpassFilter(wakeCenterZo, Fs, Fc);

% Subtract mean to get rid of previous 1D gain
wakeCenterY = wakeCenterYf - mean(wakeCenterYf);
wakeCenterZ = wakeCenterZf - mean(wakeCenterZf);
% wakeCenterY = wakeCenterYf - 0;
% wakeCenterZ = wakeCenterZf - 150;
% disp(mean(wakeCenterYf))
% disp(mean(wakeCenterZf))
% wakeCenterY = wakeCenterYf - 0;
% wakeCenterZ = wakeCenterZf - 150;

% Store wake center information in helix frame
wakecenterY_helixFrame_store = zeros(simTime, 1);
wakecenterZ_helixFrame_store = zeros(simTime, 1);
% wakecenterY_fixFrame_store = zeros(simTime, 1);
% wakecenterZ_fixFrame_store = zeros(simTime, 1);

for i = 1:1:simTime
    R_helix = [cos(omega_e*t2(i)) -sin(omega_e*t2(i)); 
               sin(omega_e*t2(i)) cos(omega_e*t2(i))];
%     invR_helix = [cos(omega_e*t2(i)) sin(omega_e*t2(i)); 
%                   -sin(omega_e*t2(i)) cos(omega_e*t2(i))];
    centerY = wakeCenterY(i);
    centerZ = wakeCenterZ(i);

    % Fixed Frame ---> Helix Frame
    center_helixFrame = R_helix * [centerZ; centerY]; % tilt, yaw
    wakecenterY_helixFrame_store(i) = center_helixFrame(2);
    wakecenterZ_helixFrame_store(i) = center_helixFrame(1);

    % Helix Frame ---> Fixed Frame
%     center_fixedFrame = invR_helix * center_helixFrame; 
%     wakecenterY_fixFrame_store(i) = center_fixedFrame(2);
%     wakecenterZ_fixFrame_store(i) = center_fixedFrame(1);
end

figure()
subplot(2, 2, 1)
plot(t, TiltYaw_FF(:, 1));
hold on;
plot(t, TiltYaw_FF(:, 2));
hold off;
title('Tilt&Yaw Fixed Frame')
legend('\beta_{tilt}', '\beta_{yaw}')

subplot(2, 2, 3)
plot(t, TiltYaw_HF(:, 1));
hold on;
plot(t, TiltYaw_HF(:, 2));
hold off;
title('Tilt&Yaw Helix Frame')
legend('\beta^e_{tilt}', '\beta^e_{yaw}')

subplot(2, 2, 2)
plot(t, wakeCenterYf)
hold on;
plot(t, wakeCenterZf)
hold off;
title('Wake Center Fixed Frame')
legend('Z - tilt', 'Y - yaw')

subplot(2, 2, 4)
plot(t, wakecenterZ_helixFrame_store)
hold on;
plot(t, wakecenterY_helixFrame_store)
hold off;
title('Wake Center Helix Frame')
legend('Z^e - tilt', 'Y^e - yaw')
end
