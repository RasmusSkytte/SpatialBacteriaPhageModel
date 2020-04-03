clearvars; close all;

path = strrep(pwd,'analysis/functions','cpp/data');

fullpath = [path '/Experiment_Conditions/beta_400_delta_0.003_eta_13200_tau_1.00_D_P_3000_nGrid_100_boundaryModifier_10_full_zeta_0.4/T_i_6.0/-1/'];

if ~exist([fullpath 'Completed.txt'],'file')
    return
end

K = [13 17];

% Load the data
Pdata = importdata([fullpath 'PhageDensity.txt']);
% Bdata = importdata([fullpath 'CellDensity.txt']);

Pdata.data = Pdata.data';
data = reshape(Pdata.data(:),100,100,4,numel(Pdata.data)/(4*100*100));

for k = 1:length(K)

% Extract the z profile just after phage addition
p{k} = squeeze(sum(sum(data(:,:,:,K(k)),1),2));

end

% Define z coordinates
z = 50:100:350;

% Create figure
fh = figure(6); clf;
fh.Position = [2640 169 2*290 290];

% Create double axes
ax1 = axes;
ax2 = axes;

% Plot the data
barh(ax1,z,p{1}/1e5,'FaceColor','r','BarWidth',1)
barh(ax2,z,p{2}/1e5,'FaceColor','r','BarWidth',1)

% Adjust the axes
title(ax1,'T = 6 h','FontSize',20)
title(ax2,'T = 8 h','FontSize',20)
ylabel(ax1,'Z (mm)','FontSize',20);
ax1.YAxis.FontSize = 20;

text(ax1,5,    -155,'P','FontSize',20,'HorizontalAlignment','center')
text(ax2,5,    -155,'P','FontSize',20,'HorizontalAlignment','center')
text(ax1,-3, 550,'b','FontSize',20,'FontWeight','bold','HorizontalAlignment','center')

ax1.YTick = [0 200 400];
ax1.TickLength = [0.05 0.01];
ax1.XLim = [0 10];
ax1.XTick = [0 5 10];
ax1.FontSize = 20;

ax1.XTickLabel = [];
text(ax1,0.8,  -150,'0','FontSize',20,'HorizontalAlignment','center')
text(ax1,9.2,  -135,'10^6','FontSize',20,'HorizontalAlignment','center')
ax1.YTickLabel = {'0.0','0.2','0.4'};

ax1.LineWidth = 2;

ax2.YTick = [];
ax2.TickLength = [0.05 0.01];
ax2.XLim = [0 10];
ax2.XTick = [0 5 10];
ax2.FontSize = 20;
ax2.XTickLabel = [];
text(ax2,0.8,  -150,'0','FontSize',20,'HorizontalAlignment','center')
text(ax2,9.2,  -135,'10^6','FontSize',20,'HorizontalAlignment','center')

ax1.Position(1) = 0.21;
ax2.Position(1) = 0.59;

ax1.Position(2) = 0.14;
ax2.Position(2) = 0.14;

ax1.Position(3) = 0.37;
ax2.Position(3) = 0.37;

ax1.Position(4) = 0.7;
ax2.Position(4) = 0.7;

ax2.LineWidth = 2;

if ~exist('../Fig_5','dir')
    mkdir('../Fig_5')
end
print('-r300','../Fig_5/upper_right.tif','-dtiff')
