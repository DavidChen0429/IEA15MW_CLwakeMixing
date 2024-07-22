function [] = MomentCenter_comparison_Visualization(SimData1, SimData2, Fs, Fc)

% Basic helix information
U_inflow = 8;        % Inflow wind speed, same with the Q-blade setting
D_IEA15MW = 240;
Str = 0.3;
Freq = Str*U_inflow/D_IEA15MW;      % From Str, in Hz
omega_e = Freq*2*pi;

LiDAR_sample1= SimData1.LiDAR_data;
TiltYaw_FF1 = SimData1.thetaTiltYaw_fixedFrame;
TiltYaw_HF1 = SimData1.thetaTiltYaw_helixFrame;
LiDAR_sample2= SimData2.LiDAR_data;
TiltYaw_FF2 = SimData2.thetaTiltYaw_fixedFrame;
TiltYaw_HF2 = SimData2.thetaTiltYaw_helixFrame;

datalength = size(LiDAR_sample1);
simTime = datalength(1); 
timeStep = 0.1;   
simLen = simTime * timeStep; % seconds
t = linspace(1, simTime, simTime);
t2 = linspace(1, simLen, simTime);

wakeCenterY1o = arrayfun(@(x) x.centerY, LiDAR_sample1);
wakeCenterZ1o = arrayfun(@(x) x.centerZ, LiDAR_sample1);
wakeCenterY2o = arrayfun(@(x) x.centerY, LiDAR_sample2);
wakeCenterZ2o = arrayfun(@(x) x.centerZ, LiDAR_sample2);

% Low pass filter
% Fs = 10;  % sampling frequency Hz
% Fc = 0.02;  % cutoff frequency Hz
wakeCenterY1f = lowpassFilter(wakeCenterY1o, Fs, Fc);
wakeCenterZ1f = lowpassFilter(wakeCenterZ1o, Fs, Fc);
wakeCenterY2f = lowpassFilter(wakeCenterY2o, Fs, Fc);
wakeCenterZ2f = lowpassFilter(wakeCenterZ2o, Fs, Fc);

% Subtract mean to get rid of previous 1D gain
wakeCenterY1 = wakeCenterY1f - mean(wakeCenterY1f);
wakeCenterZ1 = wakeCenterZ1f - mean(wakeCenterZ1o);
% Subtract mean to get rid of previous 1D gain
wakeCenterY2 = wakeCenterY2f - mean(wakeCenterY2f);
wakeCenterZ2 = wakeCenterZ2f - mean(wakeCenterZ2f);

% Store wake center information in helix frame
wakecenterY_helixFrame_store1 = zeros(simTime, 1);
wakecenterZ_helixFrame_store1 = zeros(simTime, 1);
wakecenterY_helixFrame_store2 = zeros(simTime, 1);
wakecenterZ_helixFrame_store2 = zeros(simTime, 1);

for i = 1:1:simTime
    R_helix = [cos(omega_e*t2(i)) -sin(omega_e*t2(i)); 
               sin(omega_e*t2(i)) cos(omega_e*t2(i))];
    centerY1 = wakeCenterY1(i);
    centerZ1 = wakeCenterZ1(i);
    centerY2 = wakeCenterY2(i);
    centerZ2 = wakeCenterZ2(i);

    % Fixed Frame ---> Helix Frame
    center_helixFrame1 = R_helix * [centerZ1; centerY1]; % tilt, yaw
    wakecenterY_helixFrame_store1(i) = center_helixFrame1(2);
    wakecenterZ_helixFrame_store1(i) = center_helixFrame1(1);
    center_helixFrame2 = R_helix * [centerZ2; centerY2]; % tilt, yaw
    wakecenterY_helixFrame_store2(i) = center_helixFrame2(2);
    wakecenterZ_helixFrame_store2(i) = center_helixFrame2(1);

    % Helix Frame ---> Fixed Frame
%     center_fixedFrame = invR_helix * center_helixFrame; 
%     wakecenterY_fixFrame_store(i) = center_fixedFrame(2);
%     wakecenterZ_fixFrame_store(i) = center_fixedFrame(1);
end

figure()
subplot(2, 2, 1)
plot(t, TiltYaw_FF1(:, 1), '--');
hold on;
plot(t, TiltYaw_FF1(:, 2), '--');
plot(t, TiltYaw_FF2(:, 1));
plot(t, TiltYaw_FF2(:, 2));
hold off;
title('Tilt&Yaw Fixed Frame')
legend('\beta_{tilt1}', '\beta_{yaw1}', '\beta_{tilt2}', '\beta_{yaw2}')

subplot(2, 2, 3)
plot(t, TiltYaw_HF1(:, 1), '--');
hold on;
plot(t, TiltYaw_HF1(:, 2), '--');
plot(t, TiltYaw_HF2(:, 1));
plot(t, TiltYaw_HF2(:, 2));
hold off;
title('Tilt&Yaw Helix Frame')
legend('\beta^e_{tilt1}', '\beta^e_{yaw1}','\beta^e_{tilt2}', '\beta^e_{yaw2}')

subplot(2, 2, 2)
plot(t, wakeCenterZ1f, '--')
hold on;
plot(t, wakeCenterY1f, '--')
plot(t, wakeCenterZ2f)
plot(t, wakeCenterY2f)
hold off;
title('Wake Center Fixed Frame')
legend('z_1', 'y_1', 'z_2', 'y_2')

subplot(2, 2, 4)
plot(t, wakecenterZ_helixFrame_store1, '--')
hold on;
plot(t, wakecenterY_helixFrame_store1, '--')
plot(t, wakecenterZ_helixFrame_store2)
plot(t, wakecenterY_helixFrame_store2)
hold off;
title('Wake Center Helix Frame')
legend('z_1^e', 'y_1^e', 'z_2^e', 'y_2^e')
end
