function [] = wakeCenterTraj(data1,data2,Fs,Fc)

wakeCenterY1 = arrayfun(@(x) x.centerY, data1);
wakeCenterZ1 = arrayfun(@(x) x.centerZ, data1);
wakeCenterY2 = arrayfun(@(x) x.centerY, data2);
wakeCenterZ2 = arrayfun(@(x) x.centerZ, data2);

wakeCenterY1_f = lowpassFilter(wakeCenterY1, Fs, Fc);
wakeCenterZ1_f = lowpassFilter(wakeCenterZ1, Fs, Fc);
wakeCenterY2_f = lowpassFilter(wakeCenterY2, Fs, Fc);
wakeCenterZ2_f = lowpassFilter(wakeCenterZ2, Fs, Fc);

figure('Position', [10, 10, 500, 500]);
plot(wakeCenterY1_f, wakeCenterZ1_f,'red');
hold on
plot(wakeCenterY2_f, wakeCenterZ2_f,'blue');
hold off
xlabel('Y [m]')
ylabel('Z [m]')
xlim([-50 50])
ylim([100 200])
legend('case 1', 'case 2')
title('Wake Center Trajectory')
end