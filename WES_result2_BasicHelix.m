%% Helix cases: remake 3-panel bar chart
clear; clc;

cats = categorical({'WT1','WT2','All'});
cats = reordercats(cats, {'WT1','WT2','All'});

% ---------- Replace these with your values ----------
% Power [MW]
P_S   = [3.61 1.38 4.49];
P_T   = [3.72 1.60 5.32];
P_ST  = [3.68 1.49 5.17];
P_Uni = [3.70 1.48 5.18];

% DEL Flapwise [Nm] (numbers here are in absolute units)
F_S   = [2.90e7 1.40e7 4.30e7];
F_T   = [3.00e7 1.80e7 4.80e7];
F_ST  = [3.15e7 1.82e7 4.79e7];
F_Uni = [2.80e7 1.80e7 4.60e7];

% DEL Edgewise [Nm]
E_S   = [7.50e6 3.10e6 10.60e6];
E_T   = [7.70e6 3.80e6 11.50e6];
E_ST  = [7.90e6 3.70e6 11.60e6];
E_Uni = [7.20e6 3.60e6 10.80e6];
% ---------------------------------------------------

% Colors (match MATLAB defaults used in many papers)
cS   = [0 0.4470 0.7410];     % blue
cT   = [0.8500 0.3250 0.0980];% orange
cST  = [0.3010 0.7450 0.9330];% light blue

figuresetup = [100, 100, 500, 1500];
figure('Name', 'WES Visualizaztion2', 'NumberTitle', 'off', 'Position', figuresetup);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% ------------ Panel 1: Power ------------
nexttile;
hb = bar(cats, [P_S(:) P_T(:) P_ST(:)], 'grouped', 'BarWidth', 0.84); hold on;
hb(1).FaceColor = cS;  hb(2).FaceColor = cT;  hb(3).FaceColor = cST;

% Overlay dashed-outline "Uni"
bar(cats, P_Uni(:), 'grouped', 'BarWidth', 0.84, ...
    'FaceColor','none','LineStyle','--','LineWidth',1.2,'EdgeColor','k');

ylabel('Power [MW]');
ylim([0 6]); box on; set(gca,'LineWidth',1);
legend({'Shear','Turbulence','Shear&Turbulence','Uniform'},'Location','northwest');

% ------------ Panel 2: DEL Flapwise ------------
nexttile;
scaleF = 1e7; % show as x10^7
hb = bar(cats, [F_S(:) F_T(:) F_ST(:)]/scaleF, 'grouped', 'BarWidth', 0.84); hold on;
hb(1).FaceColor = cS;  hb(2).FaceColor = cT;  hb(3).FaceColor = cST;

bar(cats, F_Uni(:)/scaleF, 'grouped', 'BarWidth', 0.84, ...
    'FaceColor','none','LineStyle','--','LineWidth',1.2,'EdgeColor','k');

ylabel('DEL Flapwise [Nm]');
text(0.03, 0.97, '\times10^7', 'Units','normalized', 'VerticalAlignment','top');
ylim([0 5.1]); box on; set(gca,'LineWidth',1);
% legend({'S','T','S\&T','Uni'},'Location','northwest');

% ------------ Panel 3: DEL Edgewise ------------
nexttile;
scaleE = 1e6; % show as x10^6
hb = bar(cats, [E_S(:) E_T(:) E_ST(:)]/scaleE, 'grouped', 'BarWidth', 0.84); hold on;
hb(1).FaceColor = cS;  hb(2).FaceColor = cT;  hb(3).FaceColor = cST;

bar(cats, E_Uni(:)/scaleE, 'grouped', 'BarWidth', 0.84, ...
    'FaceColor','none','LineStyle','--','LineWidth',1.2,'EdgeColor','k');

ylabel('DEL Edgewise [Nm]');
text(0.03, 0.97, '\times10^6', 'Units','normalized', 'VerticalAlignment','top');
ylim([0 12]); box on; set(gca,'LineWidth',1);
% legend({'S','T','S\&T','Uni'},'Location','northwest');

% Global font tweaks
set(findall(gcf,'-property','FontName'),'FontName','Times');
set(findall(gcf,'-property','FontSize'),'FontSize',11);

setfigpaper('Width',[30,0.5],'Interpreter','tex','FontSize',15,'linewidth',2)