addpath('.\Functions');

%% Beta
figure;
subplot(2, 2, 1)
plot(thetaTiltYaw_fixedFrame(:, 1));
hold on;
plot(thetaTiltYaw_fixedFrame(:, 2));
hold off;
title('\beta FF')
legend('\theta_{tilt}', '\theta_{yaw}')

[f1, P1] = FFT_func(thetaTiltYaw_fixedFrame(:, 1), 1, 10);
subplot(2, 2, 2);
plot(f1,P1,"LineWidth",1);
title("Filtered signal in Frequency Domain")
xlabel("f (Hz)")
xlim([0 0.2])
ylabel("Magnitude")
text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
annotation_text = sprintf('FFT information');
text('Units', 'normalized', 'Position', text_position, ...
    'String', annotation_text, 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'bottom');

subplot(2, 2, 3);
plot(thetaTiltYaw_helixFrame(:, 1));
hold on;
plot(thetaTiltYaw_helixFrame(: ,2));
hold off;
title('\beta_e HF')
legend('\theta^e_{tilt}', '\theta^e_{yaw}')

[f2, P2] = FFT_func(thetaTiltYaw_helixFrame(:, 1), 1, 10);
subplot(2, 2, 4)
plot(f2,P2,"LineWidth",1);
title("Filtered signal in Frequency Domain")
xlabel("f (Hz)")
xlim([0 0.2])
ylabel("Magnitude")
text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
annotation_text = sprintf('FFT information');
text('Units', 'normalized', 'Position', text_position, ...
    'String', annotation_text, 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'bottom');

%% Wake Center
t = linspace(1, simLen, simTime);
helixCenter_e = zeros(simTime, 2);
for i = 1:1:simTime
    invR_helix = [cos(omega_e*t(i)) -sin(omega_e*t(i)); 
                  sin(omega_e*t(i)) cos(omega_e*t(i))];
    centerZ = helixCenter(i, 1) - 91.9411;  % 91.9411
    centerY = helixCenter(i, 2) + 3.1245;  % -3.1245
    center_e = invR_helix * [centerZ; centerY];
    helixCenter_e(i, :) = [center_e(1) center_e(2)];
end

figure;
subplot(2, 2, 1)
plot(helixCenter(:, 1)-mean(helixCenter(:, 1)));
hold on;
plot(helixCenter(:, 2)-mean(helixCenter(:, 2)));
hold off;
title('Center FF')
legend('z', 'y')

[f11, P11] = FFT_func(helixCenter(:, 1)-mean(helixCenter(:, 1)), 1, 10);
[f12, P12] = FFT_func(helixCenter(:, 2)-mean(helixCenter(:, 2)), 1, 10);
subplot(2, 2, 2);
plot(f11,P11,"LineWidth",1);
hold on
plot(f12,P12,"LineWidth",1);
hold off
title("Filtered signal in Frequency Domain")
xlabel("f (Hz)")
xlim([0 0.2])
ylabel("Magnitude")
legend('z', 'y')
text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
annotation_text = sprintf('FFT information');
text('Units', 'normalized', 'Position', text_position, ...
    'String', annotation_text, 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'bottom');

subplot(2, 2, 3);
plot(helixCenter_e(:, 1));
hold on;
plot(helixCenter_e(: ,2));
hold off;
title('Center HF')
legend('z_e', 'y_e')

[f21, P21] = FFT_func(helixCenter_e(:, 1), 1, 10);
[f22, P22] = FFT_func(helixCenter_e(:, 2), 1, 10);
subplot(2, 2, 4);
plot(f21,P21,"LineWidth",1);
hold on
plot(f22,P22,"LineWidth",1);
hold off
title("Filtered signal in Frequency Domain")
xlabel("f (Hz)")
xlim([0 0.2])
ylabel("Magnitude")
legend('z_e', 'y_e')
text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
annotation_text = sprintf('FFT information');
text('Units', 'normalized', 'Position', text_position, ...
    'String', annotation_text, 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'bottom');