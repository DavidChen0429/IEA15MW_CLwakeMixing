%clc
%clear
close all
addpath('.\Functions');

%% Load data
Fs = 10;  % sampling frequency Hz
Fc = 0.02;  % cutoff frequency Hz
filterIndex = 500;
timeStep = 0.1;

turbineName = '.\Data\NREL5MW\';
caseName = 'HubJet\';
fileName1 = 'Baseline.mat';
fileName2 = 'St4A3.mat';
fileName3 = 'St2A3.mat';
fileName4 = 'St3A1.mat';
fileName5 = 'St3A2.mat';
fileName6 = 'St3A3.mat';
fileName7 = 'St3A4.mat';
fileName8 = 'St3A5.mat';

Baseline = load([turbineName caseName fileName1]);
St4A3 = load([turbineName caseName fileName2]);
St2A3 = load([turbineName caseName fileName3]);
St3A1 = load([turbineName caseName fileName4]);
St3A2 = load([turbineName caseName fileName5]);
St3A3 = load([turbineName caseName fileName6]);
St3A4 = load([turbineName caseName fileName7]);
St3A5 = load([turbineName caseName fileName8]);

%% Same St, Different Amplitude
t = (1: (length(Baseline.FF_beta)-filterIndex+1)) * timeStep;
lw = 2;

% Time
% z
biasZ = mean(Baseline.FF_helixCenter_filtered(filterIndex, 1)) - 90;
biasY = mean(Baseline.FF_helixCenter_filtered(filterIndex, 2)) - 0;
figure('Name', 'Time Domain', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t,Baseline.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
hold on
plot(t,St3A1.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A2.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A3.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A4.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
hold off
xlim([0 250])
ylim([70 110])
xlabel('Time [s]')
ylabel('Z [m]')
title('Z, Hub Jet Rotation with St=0.3')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4')
% setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',15,'linewidth',1)

% y
% figure('Name', 'Y', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 2)
plot(t,Baseline.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
hold on
plot(t,St3A1.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A2.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A3.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A4.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
hold off
xlim([0 250])
ylim([-20 20])
xlabel('Time [s]')
ylabel('Y [m]')
title('Y, Hub Jet Rotation with St=0.3')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4')
setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',20,'linewidth',1)

% Frequency
[f_Baseline,P1_Baseline] = FFT_func(Baseline.FF_helixCenter_filtered(:,1), filterIndex, Fs);
[f_St3A1,P1_St3A1] = FFT_func(St3A1.FF_helixCenter_filtered(:,1), filterIndex, Fs);
[f_St3A2,P1_St3A2] = FFT_func(St3A2.FF_helixCenter_filtered(:,1), filterIndex, Fs);
[f_St3A3,P1_St3A3] = FFT_func(St3A3.FF_helixCenter_filtered(:,1), filterIndex, Fs);
[f_St3A4,P1_St3A4] = FFT_func(St3A4.FF_helixCenter_filtered(:,1), filterIndex, Fs);

figure('Name', 'Frequency FFT', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
plot(f_Baseline,P1_Baseline,"LineWidth",lw);
hold on
plot(f_St3A1,P1_St3A1,"LineWidth",lw);
plot(f_St3A2,P1_St3A2,"LineWidth",lw);
plot(f_St3A3,P1_St3A3,"LineWidth",lw);
plot(f_St3A4,P1_St3A4,"LineWidth",lw);
xline(0.0238, '--k', 'LineWidth', lw);
hold off
title("Filtered signal in Frequency Domain")
xlabel("f (Hz)")
ylabel("Magnitude")
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4', 'Helix Frequency')
xlim([0 0.08])
text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
% annotation_text = sprintf('FFT information');
% text('Units', 'normalized', 'Position', text_position, ...
%     'String', annotation_text, 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'bottom');
setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',20,'linewidth',1)

%% Same Amplitude, Different St
t = (1: (length(Baseline.FF_beta)-filterIndex+1)) * timeStep;
lw = 2;

% Time
% z
biasZ = mean(Baseline.FF_helixCenter_filtered(filterIndex, 1)) - 90;
biasY = mean(Baseline.FF_helixCenter_filtered(filterIndex, 2)) - 0;
figure('Name', 'Time Domain', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t,Baseline.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
hold on
plot(t,St2A3.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A3.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St4A3.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
hold off
xlim([0 250])
ylim([70 110])
xlabel('Time [s]')
ylabel('Z [m]')
title('Z, Hub Jet Rotation with A=3')
legend('Baseline', 'St=0.2', 'St=0.3', 'St=0.4')
% setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',15,'linewidth',1)

% y
% figure('Name', 'Y', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 2)
plot(t,Baseline.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
hold on
plot(t,St2A3.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A3.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St4A3.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
hold off
xlim([0 250])
ylim([-20 20])
xlabel('Time [s]')
ylabel('Y [m]')
title('Y, Hub Jet Rotation with A=3')
legend('Baseline', 'St=0.2', 'St=0.3', 'St=0.4')
setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',20,'linewidth',1)

% Frequency
[f_Baseline,P1_Baseline] = FFT_func(Baseline.FF_helixCenter_filtered(:,1), filterIndex, Fs);
[f_St4A3,P1_St4A3] = FFT_func(St4A3.FF_helixCenter_filtered(:,1), filterIndex, Fs);
[f_St2A3,P1_St2A3] = FFT_func(St2A3.FF_helixCenter_filtered(:,1), filterIndex, Fs);
[f_St3A3,P1_St3A3] = FFT_func(St3A3.FF_helixCenter_filtered(:,1), filterIndex, Fs);

figure('Name', 'Frequency FFT', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
plot(f_Baseline,P1_Baseline,"LineWidth",lw);
hold on
plot(f_St2A3,P1_St2A3,"LineWidth",lw);
plot(f_St3A3,P1_St3A3,"LineWidth",lw);
plot(f_St4A3,P1_St4A3,"LineWidth",lw);
a = xline(0.2*10/126, '--k', 'LineWidth', lw);
a.Color = "#A2142F";
c = xline(0.3*10/126, '--k', 'LineWidth', lw);
c.Color = "#EDB120";
d = xline(0.4*10/126, '--r', 'LineWidth', lw);
d.Color = "#7E2F8E";
hold off
title("Filtered signal in Frequency Domain")
xlabel("f (Hz)")
ylabel("Magnitude")
legend('Baseline', 'St=0.2', 'St=0.3', 'St=0.4', 'Helix1','Helix2','Helix3')
xlim([0 0.08])
ylim([-0.5 7])
% text_position = [0.95, 0.05]; % Normalized coordinates (bottom right)
% annotation_text = sprintf('FFT information');
% text('Units', 'normalized', 'Position', text_position, ...
%     'String', annotation_text, 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'bottom');
setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',20,'linewidth',1)

%% Helix Frame Transform
t = (1: (length(Baseline.FF_beta)-filterIndex+1)) * timeStep;
lw = 2;

% === Output
% Fixed Frame
% z
biasZ = mean(Baseline.FF_helixCenter_filtered(filterIndex, 1)) - 90;
biasY = mean(Baseline.FF_helixCenter_filtered(filterIndex, 2)) - 0;
figure('Name', 'Fixed Frame', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t,Baseline.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
hold on
plot(t,St3A1.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A2.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A3.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A4.FF_helixCenter_filtered(filterIndex:end, 1)-biasZ,'LineWidth', lw)
hold off
xlim([0 250])
ylim([70 110])
xlabel('Time [s]')
ylabel('Z [m]')
title('Z  Fixed Frame')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4', 'Location', 'southeast')
% setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',15,'linewidth',1)
% y
% figure('Name', 'Y', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 2)
plot(t,Baseline.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
hold on
plot(t,St3A1.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A2.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A3.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A4.FF_helixCenter_filtered(filterIndex:end, 2)-biasY,'LineWidth', lw)
hold off
xlim([0 250])
ylim([-20 20])
xlabel('Time [s]')
ylabel('Y [m]')
title('Y  Fixed Frame')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4', 'Location', 'southeast')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',1)

% Helix Frame
figure('Name', 'Helix Frame', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t,Baseline.HF_helixCenter_filtered(filterIndex:end, 1),'LineWidth', lw)
hold on
plot(t,St3A1.HF_helixCenter_filtered(filterIndex:end, 1),'LineWidth', lw)
plot(t,St3A2.HF_helixCenter_filtered(filterIndex:end, 1),'LineWidth', lw)
plot(t,St3A3.HF_helixCenter_filtered(filterIndex:end, 1),'LineWidth', lw)
plot(t,St3A4.HF_helixCenter_filtered(filterIndex:end, 1),'LineWidth', lw)
hold off
xlim([0 250])
xlabel('Time [s]')
ylabel('Z [m]')
title('Z^e  Helix Frame')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4', 'Location', 'southeast')
% setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',15,'linewidth',1)
% y
% figure('Name', 'Y', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 2)
plot(t,Baseline.HF_helixCenter_filtered(filterIndex:end, 2),'LineWidth', lw)
hold on
plot(t,St3A1.HF_helixCenter_filtered(filterIndex:end, 2),'LineWidth', lw)
plot(t,St3A2.HF_helixCenter_filtered(filterIndex:end, 2),'LineWidth', lw)
plot(t,St3A3.HF_helixCenter_filtered(filterIndex:end, 2),'LineWidth', lw)
plot(t,St3A4.HF_helixCenter_filtered(filterIndex:end, 2),'LineWidth', lw)
hold off
xlim([0 250])
xlabel('Time [s]')
ylabel('Y [m]')
title('Y^e  Helix Frame')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4', 'Location', 'southeast')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',1)


% === Input
% Fixed Frame
% z
biasZ = 0;
biasY = 0;
figure('Name', 'Fixed Frame', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t,Baseline.FF_beta(filterIndex:end, 1)-biasZ,'LineWidth', lw)
hold on
plot(t,St3A1.FF_beta(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A2.FF_beta(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A3.FF_beta(filterIndex:end, 1)-biasZ,'LineWidth', lw)
plot(t,St3A4.FF_beta(filterIndex:end, 1)-biasZ,'LineWidth', lw)
hold off
xlim([0 250])
% ylim([70 110])
xlabel('Time [s]')
ylabel('Magnitude')
title('\beta_{tilt}  Fixed Frame')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4', 'Location', 'southeast')
% setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',15,'linewidth',1)
% y
% figure('Name', 'Y', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 2)
plot(t,Baseline.FF_beta(filterIndex:end, 2)-biasY,'LineWidth', lw)
hold on
plot(t,St3A1.FF_beta(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A2.FF_beta(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A3.FF_beta(filterIndex:end, 2)-biasY,'LineWidth', lw)
plot(t,St3A4.FF_beta(filterIndex:end, 2)-biasY,'LineWidth', lw)
hold off
xlim([0 250])
% ylim([-20 20])
xlabel('Time [s]')
ylabel('Magnitude')
title('\beta_{yaw}  Fixed Frame')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4', 'Location', 'southeast')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',1)

% Helix Frame
figure('Name', 'Helix Frame', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 1)
plot(t,Baseline.HF_beta(filterIndex:end, 1),'LineWidth', lw)
hold on
plot(t,St3A1.HF_beta(filterIndex:end, 1),'LineWidth', lw)
plot(t,St3A2.HF_beta(filterIndex:end, 1),'LineWidth', lw)
plot(t,St3A3.HF_beta(filterIndex:end, 1),'LineWidth', lw)
plot(t,St3A4.HF_beta(filterIndex:end, 1),'LineWidth', lw)
hold off
xlim([0 250])
xlabel('Time [s]')
ylabel('Magnitude')
title('\beta^e_{tilt}  Helix Frame')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4', 'Location', 'southeast')
% setfigpaper('Width',[30,0.5],'Interpreter','Latex','FontSize',15,'linewidth',1)
% y
% figure('Name', 'Y', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 600]);
subplot(2, 1, 2)
plot(t,Baseline.HF_beta(filterIndex:end, 2),'LineWidth', lw)
hold on
plot(t,St3A1.HF_beta(filterIndex:end, 2),'LineWidth', lw)
plot(t,St3A2.HF_beta(filterIndex:end, 2),'LineWidth', lw)
plot(t,St3A3.HF_beta(filterIndex:end, 2),'LineWidth', lw)
plot(t,St3A4.HF_beta(filterIndex:end, 2),'LineWidth', lw)
hold off
xlim([0 250])
xlabel('Time [s]')
ylabel('Magnitude')
title('\beta^e_{yaw}  Helix Frame')
legend('Baseline', 'A=1', 'A=2', 'A=3', 'A=4', 'Location', 'southeast')
setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',20,'linewidth',1)