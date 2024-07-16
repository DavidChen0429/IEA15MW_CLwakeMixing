%% This file is created to study the correlation between tilt,yaw to wake center
clear
close all
addpath('.\Functions');

%% Basic information definition
fileName = '600s_Center_HF_Mtilte.mat';
dataPath = '.\Data\MAT\LiDAR_sampling\';
caseName = 'Uni\Str0.3_U8_1Dd_10Hz_CCW\';
SimData = load([dataPath caseName fileName]);
wakeCenterChange_Visualization(SimData)
ringVisualization(SimData.LiDAR_data)

% LiDAR_sample= SimData.LiDAR_data;
% TiltYaw_FF = SimData.thetaTiltYaw_fixedFrame;
% TiltYaw_HF = SimData.thetaTiltYaw_helixFrame;
% 
% datalength = size(LiDAR_sample);
% t = linspace(1, datalength(2), datalength(2));
% wakeCenterY = arrayfun(@(x) x.centerY, LiDAR_sample);
% wakeCenterZ = arrayfun(@(x) x.centerZ, LiDAR_sample);
% 
% figure()
% subplot(3, 1, 1)
% plot(t, TiltYaw_FF(:, 1));
% hold on;
% plot(t, TiltYaw_FF(:, 2));
% hold off;
% title('Fixed Frame')
% legend('\beta_{tilt}', '\beta_{yaw}')
% 
% subplot(3, 1, 2)
% plot(t, TiltYaw_HF(:, 1));
% hold on;
% plot(t, TiltYaw_HF(:, 2));
% hold off;
% title('Helix Frame')
% legend('\beta^e_{tilt}', '\beta^e_{yaw}')
% 
% subplot(3, 1, 3)
% plot(t, wakeCenterZ)
% hold on;
% plot(t, wakeCenterY)
% hold off;
% title('Wake Center')
% legend('Z - tilt', 'Y - yaw')