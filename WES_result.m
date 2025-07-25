clear
close all

% Example dummy data
% Dimensions: [controlType, turbine, subplotIndex]
% controlType = 1: Open-loop (OL), 2: Closed-loop (CL)
data = rand(2, 3, 9);  % [OL/CL x WT1/WT2/All x 9 subplots]

% 1. Change value of data
% 2. Fix scale

% Assigning value to data array
% ======== Shear
% subplot(1,1) --- power
data(:, :, 1) = [3.65, 1.4, 5.05;   % OL (WT1,WT2,All)
                 3.6, 1.46, 5.06];  % CL (WT1,WT2,All)
% subplot(2,1) --- DELf
data(:, :, 4) = [2.9e7, 1.45e7, 4.35e7;   % OL (WT1,WT2,All)
                 3.15e7, 1.75e7, 4.9e7];  % CL (WT1,WT2,All)
% subplot(3,1) --- DELe
data(:, :, 7) = [7.4e6, 3.2e6, 10.6e6;   % OL (WT1,WT2,All)
                 8.3e6, 3.8e6, 12.1e6];  % CL (WT1,WT2,All)
% ======== Turbulence
% subplot(1,2) --- power
data(:, :, 2) = [3.70, 1.55, 5.25;   % OL (WT1,WT2,All)
                 3.71, 1.55, 5.26];  % CL (WT1,WT2,All)
% subplot(2,2) --- DELf
data(:, :, 5) = [3e7, 1.85e7, 4.85e7;   % OL (WT1,WT2,All)
                 2.95e7, 1.95e7, 4.9e7];  % CL (WT1,WT2,All)
% subplot(3,2) --- DELe
data(:, :, 8) = [7.6e6, 3.7e6, 11.3e6;   % OL (WT1,WT2,All)
                 7.59e6, 3.71e6, 11.3e6];  % CL (WT1,WT2,All)

% ======== Combined
% subplot(1,3) --- power
data(:, :, 3) = [3.68, 1.51, 5.19;   % OL (WT1,WT2,All)
                 3.51, 1.70, 5.21];  % CL (WT1,WT2,All)
% subplot(2,3) --- DELf
data(:, :, 6) = [3.2e7, 1.87e7, 5.07e7;   % OL (WT1,WT2,All)
                 3.9e7, 1.90e7, 5.80e7];  % CL (WT1,WT2,All)
% subplot(3,3) --- DELe
data(:, :, 9) = [7.9e6, 3.60e6, 11.50e6;   % OL (WT1,WT2,All)
                 8.5e6, 3.96e6, 12.46e6];  % CL (WT1,WT2,All)

blue  = [0, 0.4470, 0.7410];    % OL
orange = [0.8500, 0.3250, 0.0980];  % CL

% Labels
colTitles = {'Shear', 'Turbulence', 'Combined'};
rowLabels = {'Power [MW]', 'DEL Flapwise [Nm]', 'DEL Edgewise [Nm]'};
xLabels = {'WT1', 'WT2', 'All'};

figuresetup = [100, 100, 500, 400];
figure('Name', 'WES Visualizaztion', 'NumberTitle', 'off', 'Position', figuresetup);
t = tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:9
    ax = nexttile(i);
    hold on;

    % Extract values for this subplot
    y_OL = data(1, :, i);  % Open-loop 
    y_CL = data(2, :, i);  % Closed-loop 

    % Bar parameters
    x = 1:3;
    bw = 0.35;

    % Plot Open-loop bars 
    bar(x - bw/2, y_OL, bw, 'FaceColor', blue, 'EdgeColor', 'k');
    % Plot Closed-loop bars 
    bar(x + bw/2, y_CL, bw, 'FaceColor', orange, 'EdgeColor', 'k');

    % Add annotations: percentage increase from OL to CL
    for j = 1:3  % Loop over WT1, WT2, All
        pct = (y_CL(j) - y_OL(j)) / y_OL(j) * 100;
        text_x = x(j);  % Center between bars
        if j == 3
            text_y = max(y_CL(j), y_OL(j)) * 1.1;  % Slightly above taller bar
        else
            text_y = max(y_CL(j), y_OL(j)) * 1.2;  % Slightly above taller bar
        end
        % Format string: show + or - with one decimal
        pctStr = sprintf('%+.1f%%', pct);

        % Change text color
        if pct > 0
            txtColor = 'green';
        else
            txtColor = 'red';
        end

        % Annotate
        text(text_x, text_y, pctStr, ...
             'Color', txtColor, ... 
             'HorizontalAlignment', 'center', ...
             'FontSize', 9, ...
             'FontWeight', 'bold');
    end

    % Formatting
    xlim([0.5 3.5]);
    ylim([0, max([y_OL y_CL]) * 1.2]);
    if i > 6
        xticks(x);
        xticklabels(xLabels);
    else
        xticks([]);       % Hide ticks
        xticklabels({});  % Hide labels
    end
    grid on;

    % Add column titles to top row
    if i <= 3
        title(colTitles{i});
    end

    % Add row labels to first column
    if mod(i-1, 3) == 0
        ylabel(rowLabels{(i-1)/3 + 1});
    end

    % Legend only in first subplot
    if i == 1
        h1 = bar(nan, nan, 'FaceColor', blue, 'EdgeColor', 'k');
        h2 = bar(nan, nan, 'FaceColor', orange, 'EdgeColor', 'k');
        legend([h1 h2], {'OL', 'CL'}, 'Location', 'southwest');
    end

    hold off;
end
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',15,'linewidth',2)
    