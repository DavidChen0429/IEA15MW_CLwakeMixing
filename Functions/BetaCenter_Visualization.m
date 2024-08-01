function [] = BetaCenter_Visualization(SimData, Fs, Fc)
len = length(SimData.LiDAR_data);
t = linspace(1, len*0.1, len);

figure()
subplot(2, 2, 1)
plot(t, SimData.FF_theta(:, 1));
hold on;
plot(t, SimData.FF_theta(:, 2));
hold off;
title('Tilt&Yaw Fixed Frame')
legend('\beta_{tilt}', '\beta_{yaw}')

subplot(2, 2, 3)
plot(t, SimData.HF_theta(:, 1));
hold on;
plot(t, SimData.HF_theta(:, 2));
hold off;
title('Tilt&Yaw Helix Frame')
legend('\beta^e_{tilt}', '\beta^e_{yaw}')

subplot(2, 2, 2)
plot(t, SimData.FF_helixCenter(:, 1))
hold on;
plot(t, SimData.FF_helixCenter(:, 2))
hold off;
title('Wake Center Fixed Frame')
legend('z', 'y')

subplot(2, 2, 4)
plot(t, lowpassFilter(SimData.HF_helixCenter(:, 1), Fs, Fc))
hold on;
plot(t, lowpassFilter(SimData.HF_helixCenter(:, 2), Fs, Fc))
hold off;
title('Wake Center Helix Frame')
legend('z_e', 'y_e')
end
