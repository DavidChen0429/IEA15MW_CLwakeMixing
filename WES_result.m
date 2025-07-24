% Example dummy data
% Dimensions: [controlType, turbine, subplotIndex]
% controlType = 1: Open-loop (OL), 2: Closed-loop (CL)
data = rand(2, 3, 9);  % [OL/CL x WT1/WT2/All x 9 subplots]

% Assigning value to data array
% ======== Shear
% subplot(1,1) --- shear&power
data(:, :, 1) = [1.2, 0.9, 2.8;   % OL (WT1,WT2,All)
                 1.5, 1.1, 3.2];  % CL (WT1,WT2,All)
% subplot(2,1) --- shear&DELf
data(:, :, 4) = [1e5, 0.9e5, 1.4e5;   % OL (WT1,WT2,All)
                 1.5e5, 1.1e5, 3.2e5];  % CL (WT1,WT2,All)
% subplot(3,1) --- shear&DELe
data(:, :, 7) = [1e5, 0.9e5, 1.4e5;   % OL (WT1,WT2,All)
                 1.5e5, 1.1e5, 3.2e5];  % CL (WT1,WT2,All)
% ======== Turbulence
% ======== Combined

blue  = [0, 0.4470, 0.7410];    % OL
orange = [0.8500, 0.3250, 0.0980];  % CL

% Labels
colTitles = {'Shear', 'Turbulence', 'Combined'};
rowLabels = {'Power [MW]', 'DEL Flapwise [Nm]', 'DEL Edgewise [Nm]'};
xLabels = {'WT1', 'WT2', 'All'};

% figure('Position', [100, 100, 500, 400]);
figure('Name', 'WES Visualizaztion', 'NumberTitle', 'off', 'Position', [100, 100, 500, 400]);
t = tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:9
    ax = nexttile(i);
    hold on;

    % Extract values for this subplot
    y_OL = data(1, :, i);  % Open-loop (red)
    y_CL = data(2, :, i);  % Closed-loop (green)

    % Bar parameters
    x = 1:3;
    bw = 0.35;

    % Plot Open-loop bars (red)
    bar(x - bw/2, y_OL, bw, 'FaceColor', blue, 'EdgeColor', 'k');

    % Plot Closed-loop bars (green)
    bar(x + bw/2, y_CL, bw, 'FaceColor', orange, 'EdgeColor', 'k');

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
        legend([h1 h2], {'OL', 'CL'}, 'Location', 'northwest');
    end

    hold off;
end
    setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',15,'linewidth',2)
    