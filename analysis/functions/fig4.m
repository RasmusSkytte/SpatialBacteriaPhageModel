clearvars; close all;

path = strrep(pwd,'analysis/functions','cpp/data');

mar = 50;
sep = 25;
hei = 200;
wid = 200;

% number of simulations
n = 4;

% Prepare figure
fh = figure();
fh.Position = [750 50 3*wid+2*mar+3*sep+2*sep n*hei+2*mar+n*sep];


clf
% Prepare axes;
ax = gobjects(n,3);
ch = gobjects(n,3);

for i = 1:n

    subpaths = dir(sprintf('%s/Model_%d/',path,i));
    % Load data
    switch i
        case 1
            importPath = sprintf('%s/Model_%d/%s/data_Model_%d.mat',path,i,subpaths(end).name,i);
        case 2
            importPath = sprintf('%s/Model_%d/%s/data_Model_%d.mat',path,i,subpaths(end).name,i);
        case 3
            importPath = sprintf('%s/Model_%d/%s/data_Model_%d.mat',path,i,subpaths(end).name,i);
        case 4
            importPath = sprintf('%s/Model_%d/%s/data_Model_%d.mat',path,i,subpaths(11).name,i);
    end
    model = importdata(importPath);


    T = [5 10 15];
    for j = 1:3

        % Compute index
        I = 3*(i-1)+j;

        % Create new axes
        ax(I) = axes;
        ax(I).Units = 'pixels';
        ax(I).XTick = [];
        ax(I).YTick = [];

        ax(I).OuterPosition = [mar+sep+(j-1)*(wid+sep)    (4-i)*(hei+sep)+mar+sep   wid hei];
        ax(I).Position      = [mar+sep+(j-1)*(wid+sep)    (4-i)*(hei+sep)+mar+sep   wid hei];

        % Calculate which time to plot:
        t = T(j);

        % Plot data
        ind = model.T == t;

        PP = log10(model.P(:,:,ind));
        BB = log10(model.B(:,:,ind));
        II = log10(model.I(:,:,ind));

        imagesc(ax(I),log10(model.b(:,:,ind))); hold on;
        contour(ax(I),log10(model.n(:,:,ind)),log10((1e9-2e8)*[1 1]),'r','LineWidth',2);
        ch(I) = ax(I).Children(1);
        scatter(ax(I),17.6,11,'wo','filled');

        % Set axes options
        ax(I).YDir = 'normal';
        ax(I).CLim = [0 9];
        ax(I).XLim = [1-0.5 sum(model.PFU <= 1e9)+0.5];
        ax(I).FontSize = 14;

        ax(I).LineWidth = 2;

        if j > 1
            ax(I).YTick = [];
        else
            ax(I).YTick = 1:6:ax(I).YLim(2);
            ax(I).YTickLabel = cellfun(@(x)sprintf('10^{%d}',round(log10(model.CFU(x)),2)),num2cell(ax(I).YTick),'UniformOutput',false);
        end
        if i < 4
            ax(I).XTick = [];
        else
            ax(I).XTick = 1:6:ax(I).XLim(2);
            ax(I).XTickLabel = cellfun(@(x)sprintf('10^{%d}',round(log10(model.PFU(x)),2)),num2cell(ax(I).XTick),'UniformOutput',false);
            ax(I).XTick(1)   = ax(I).XTick(1)+0.5;
            ax(I).XTick(end) = ax(I).XTick(end)-0.5;
        end

        if i == 1
            title(sprintf('T = %d h',j*5),'FontWeight','Normal')
        end
    end
end

% Add labels to graph
ax0 = axes('Position',[0 0 1 1],'Visible','off');

% X Label
hx = text(0.5,0.015,'Initial phage density (PFU / ml)');
hx.FontSize = 14;
% hx.FontWeight = 'bold';
hx.HorizontalAlignment = 'center';

% Y Label
hy = text(0.02,0.5,'Initial bacterial density (CFU / ml)');
hy.Rotation = 90;
hy.FontSize = 14;
% hy.FontWeight = 'bold';
hy.HorizontalAlignment = 'center';

% Cbar Label
hcbar = text(0.97,0.515,'Bacterial density at time T (CFU / ml)');
hcbar.Rotation = 90;
hcbar.FontSize = 14;
% hy.FontWeight = 'bold';
hcbar.HorizontalAlignment = 'center';

% A/B/C/D Label
hA = text(mar+1.5*sep,3*(hei+sep)+hei+mar,'a','Units','Pixels');
hA.Color = 'w';
hA.FontWeight = 'bold';
hA.FontSize = 16;

hB = text(mar+1.5*sep,2*(hei+sep)+hei+mar,'b','Units','Pixels');
hB.Color = 'k';
hB.FontWeight = 'bold';
hB.FontSize = 16;

hC = text(mar+1.5*sep, (hei+sep)+hei+mar,'c','Units','Pixels');
hC.Color = 'k';
hC.FontWeight = 'bold';
hC.FontSize = 16;

hD = text(mar+1.5*sep,  hei+mar,'d','Units','Pixels');
hD.Color = 'k';
hD.FontWeight = 'bold';
hD.FontSize = 16;

% Add colorbar
cb = colorbar(ax0);
% cb.Location = 'south';
% cb.Position = [m 0.95 1-m-w 0.015];
cb.Units = 'pixels';
cb.Position = [mar+3*sep+3*wid+sep/2 mar+sep 20 hei*4+3*sep];
cb.Limits = [0 10];
cb.FontSize = 12;
cb.Ticks = 0:2:10;
cb.TickLabels = cellfun(@(x)sprintf('10^{%d}',x),num2cell(cb.Ticks'),'UniformOutput',false);
caxis([0 10])



delete(ch(1:3))


fh.Color = [1 1 1];
set(gcf,'PaperPositionMode','auto')
set(gcf,'InvertHardcopy','off')

if ~exist('../Fig_4','dir')
    mkdir('../Fig_4')
end
print('-r300','../Fig_4/Fig4.tif','-dtiff')


