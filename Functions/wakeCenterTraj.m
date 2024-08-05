function [] = wakeCenterTraj(data1,data2,Fs,Fc)

wakeCenterZ1 = data1.FF_helixCenter(:, 1);
wakeCenterY1 = data1.FF_helixCenter(:, 2);
wakeCenterZ2 = data2.FF_helixCenter(:, 1);
wakeCenterY2 = data2.FF_helixCenter(:, 2);

wakeCenterY1_f = lowpassFilter(wakeCenterY1, Fs, Fc);
wakeCenterZ1_f = lowpassFilter(wakeCenterZ1, Fs, Fc);
wakeCenterY2_f = lowpassFilter(wakeCenterY2, Fs, Fc);
wakeCenterZ2_f = lowpassFilter(wakeCenterZ2, Fs, Fc);

figure('Position', [10, 10, 500, 500]);
plot(wakeCenterY1_f, wakeCenterZ1_f,'red');
hold on
% plot(mean(wakeCenterY1_f), mean(wakeCenterZ1_f),'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(wakeCenterY2_f, wakeCenterZ2_f,'blue');
% plot(mean(wakeCenterY2_f), mean(wakeCenterZ2_f),'bo', 'MarkerSize', 10, 'LineWidth', 2);
hold off
xlabel('Y [m]')
ylabel('Z [m]')
% text(mean(wakeCenterY1_f), mean(wakeCenterZ1_f), sprintf('(%0.2f, %0.2f)', mean(wakeCenterY1_f), mean(wakeCenterZ1_f)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
% text(mean(wakeCenterY2_f), mean(wakeCenterZ2_f), sprintf('(%0.2f, %0.2f)', mean(wakeCenterY2_f), mean(wakeCenterZ2_f)), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
xlim([-30 30])
ylim([60 120])
legend('Traj1','Traj2')
% legend('Traj1', 'Cetner1','Traj2', 'Cetner2')
title('Wake Center Trajectory')
end