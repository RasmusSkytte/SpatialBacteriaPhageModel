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
ax = [];
for i = 1:n
    
    % Load data
    switch i
        case 1
            subpaths = dir(sprintf('%s/Model_3',path));
            importPath = sprintf('%s/Model_3/%s/data_Model_3.mat',path,subpaths(end).name);
        case 2
            subpaths = dir(sprintf('%s/Model_4_onlyClustering',path));
            importPath = sprintf('%s/Model_4_onlyClustering/%s/data_Model_4_onlyClustering.mat',path,subpaths(end).name);
        case 3
            subpaths = dir(sprintf('%s/Model_4_onlyShielding',path));
            importPath = sprintf('%s/Model_4_onlyShielding/%s/data_Model_4_onlyShielding.mat',path,subpaths(end).name);
        case 4
            subpaths = dir(sprintf('%s/Model_4_onlyReadsorption',path));
            importPath = sprintf('%s/Model_4_onlyReadsorption/%s/data_Model_4_onlyReadsorption.mat',path,subpaths(end).name);
    end
    model = importdata(importPath);
        
    T = [5 10 15];
    for j = 1:3
        
        % Compute index
        I = 3*(i-1)+j;
        
        % Create new axes
        ax = [ax axes];
        ax(I).Units = 'pixels';
        ax(I).XTick = [];
        ax(I).YTick = [];
        
        ax(I).OuterPosition = [mar+sep+(j-1)*(wid+sep)    (n-i)*(hei+sep)+mar+sep   wid hei];
        ax(I).Position      = [mar+sep+(j-1)*(wid+sep)    (n-i)*(hei+sep)+mar+sep   wid hei];
        
        % Calculate which time to plot:
        t = T(j);
        
        % Plot data
        if size(model.b,1) ~= numel(model.CFU)
            imagesc(ax(I),log10(model.b(:,:,model.T == t)'));
        else
            imagesc(ax(I),log10(model.b(:,:,model.T == t)));
        end
        
        
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
        if i < n
            ax(I).XTick = [];
        else
            ax(I).XTick = 1:6:ax(I).XLim(2);
            ax(I).XTickLabel = cellfun(@(x)sprintf('10^{%d}',round(log10(model.PFU(x)),2)),num2cell(ax(I).XTick),'UniformOutput',false);
            ax(I).XTick(1)   = ax(I).XTick(1)+0.5;
            ax(I).XTick(end) = ax(I).XTick(end)-0.5;
        end
        
        if i == 1
            title(sprintf('T = %d h',j*5),'FontWeight','Normal','FontSize',22)
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
hA = text(mar+1.5*sep,(n-1)*(hei+sep)+hei+mar,'c','Units','Pixels');
hA.Color = 'k';
hA.FontWeight = 'bold';
hA.FontSize = 16;

hA = text(mar+1.5*sep,(n-2)*(hei+sep)+hei+mar,'c','Units','Pixels');
hA.Color = 'k';
hA.FontWeight = 'bold';
hA.FontSize = 16;

hA = text(mar+1.2*sep,(n-2)*(hei+sep)-sep+hei+mar,'+ Clustering','Units','Pixels');
hA.Color = 'k';
hA.FontWeight = 'bold';
hA.FontSize = 10;

hB = text(mar+1.5*sep,(n-3)*(hei+sep)+hei+mar,'c','Units','Pixels');
hB.Color = 'k';
hB.FontWeight = 'bold';
hB.FontSize = 16;

hB = text(mar+1.2*sep,(n-3)*(hei+sep)-sep+hei+mar,'+ Shielding','Units','Pixels');
hB.Color = 'k';
hB.FontWeight = 'bold';
hB.FontSize = 10;

hC = text(mar+1.5*sep, (n-4)*(hei+sep)+hei+mar,'c','Units','Pixels');
hC.Color = 'k';
hC.FontWeight = 'bold';
hC.FontSize = 16;

hC = text(mar+1.2*sep, (n-4)*(hei+sep)-sep+hei+mar,'+ Readsorb','Units','Pixels');
hC.Color = 'k';
hC.FontWeight = 'bold';
hC.FontSize = 10;

% Add colorbar
cb = colorbar(ax0);
% cb.Location = 'south';
% cb.Position = [m 0.95 1-m-w 0.015];
cb.Units = 'pixels';
cb.Position = [mar+3*sep+3*wid+sep/2 mar+sep 20 hei*n+(n-1)*sep];
cb.Limits = [0 10];
cb.FontSize = 12;
cb.Ticks = 0:2:10;
cb.TickLabels = cellfun(@(x)sprintf('10^{%d}',x),num2cell(cb.Ticks'),'UniformOutput',false);
caxis([0 10])


fh.Color = [1 1 1];
set(gcf,'PaperPositionMode','auto')
set(gcf,'InvertHardcopy','off')


if ~exist('../Fig_S6','dir')
    mkdir('../Fig_S6')
end
print('../Fig_S6/FigS6.tif','-dtiff')