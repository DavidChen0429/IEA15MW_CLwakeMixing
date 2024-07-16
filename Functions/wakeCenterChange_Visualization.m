function [] = wakeCenterChange_Visualization(SimData)
LiDAR_sample= SimData.LiDAR_data;
TiltYaw_FF = SimData.thetaTiltYaw_fixedFrame;
TiltYaw_HF = SimData.thetaTiltYaw_helixFrame;

datalength = size(LiDAR_sample);
t = linspace(1, datalength(1), datalength(1));
wakeCenterY = arrayfun(@(x) x.centerY, LiDAR_sample);
wakeCenterZ = arrayfun(@(x) x.centerZ, LiDAR_sample);

figure()
subplot(3, 1, 1)
plot(t, TiltYaw_FF(:, 1));
hold on;
plot(t, TiltYaw_FF(:, 2));
hold off;
title('Fixed Frame')
legend('\beta_{tilt}', '\beta_{yaw}')

subplot(3, 1, 2)
plot(t, TiltYaw_HF(:, 1));
hold on;
plot(t, TiltYaw_HF(:, 2));
hold off;
title('Helix Frame')
legend('\beta^e_{tilt}', '\beta^e_{yaw}')

subplot(3, 1, 3)
plot(t, wakeCenterZ)
hold on;
plot(t, wakeCenterY)
hold off;
title('Wake Center')
legend('Z - tilt', 'Y - yaw')
end
