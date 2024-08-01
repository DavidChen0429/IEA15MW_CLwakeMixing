function [] = BetaCenter_Comparison_Visualization(SimData1, SimData2, Fs, Fc)
len = length(SimData1.LiDAR_data);
t = linspace(1, len*0.1, len);

figure()
subplot(2, 2, 1)
plot(t, SimData1.FF_theta(:, 1), '--');
hold on;
plot(t, SimData1.FF_theta(:, 2), '--');
plot(t, SimData2.FF_theta(:, 1));
plot(t, SimData2.FF_theta(:, 2));
hold off;
title('Tilt&Yaw Fixed Frame')
legend('\beta_{tilt1}', '\beta_{yaw1}', '\beta_{tilt2}', '\beta_{yaw2}')

subplot(2, 2, 3)
plot(t, SimData1.HF_theta(:, 1), '--');
hold on;
plot(t, SimData1.HF_theta(:, 2), '--');
plot(t, SimData2.HF_theta(:, 1));
plot(t, SimData2.HF_theta(:, 2));
hold off;
title('Tilt&Yaw Helix Frame')
legend('\beta^e_{tilt1}', '\beta^e_{yaw1}','\beta^e_{tilt2}', '\beta^e_{yaw2}')

subplot(2, 2, 2)
plot(t, SimData1.FF_helixCenter(:, 1), '--')
hold on;
plot(t, SimData1.FF_helixCenter(:, 2), '--')
plot(t, SimData2.FF_helixCenter(:, 1))
plot(t, SimData2.FF_helixCenter(:, 2))
hold off;
title('Wake Center Fixed Frame')
legend('z_1', 'y_1', 'z_2', 'y_2')

subplot(2, 2, 4)
plot(t, lowpassFilter(SimData1.HF_helixCenter(:, 1), Fs, Fc), '--')
hold on;
plot(t, lowpassFilter(SimData1.HF_helixCenter(:, 2), Fs, Fc), '--')
plot(t, lowpassFilter(SimData2.HF_helixCenter(:, 1), Fs, Fc))
plot(t, lowpassFilter(SimData2.HF_helixCenter(:, 2), Fs, Fc))
hold off;
title('Wake Center Helix Frame')
legend('z_1^e', 'y_1^e', 'z_2^e', 'y_2^e')
end
