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


nGrid = [10 25 50];
dT = [2e-3 1e-3 5e-4 2e-4];

clf
% Prepare axes;
ax = [];
for i = 1:n
    
    T = 15;
    for j = 1:3

        % Load data
        importPath = sprintf('%s/Model_4/beta_100_delta_0.100_eta_10000_tau_0.50', path);
        if i == 1
            importPath = sprintf('%s_nGrid_%d/data_Model_4.mat', importPath, nGrid(j));
        else
            importPath = sprintf('%s_nGrid_%d_dt_%.2fx/data_Model_4.mat', importPath, nGrid(j), dT(i)/2e-3);
        end
        model = importdata(importPath);

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
        t = T;
        
        % Truncate data
        model.b = model.b(model.PFU >= 1e3,model.CFU >= 1e3,:);
        model.PFU = model.PFU(model.PFU >= 1e3);
        model.CFU = model.CFU(model.CFU >= 1e3);
        
        % Plot data
        imagesc(ax(I),log10(model.b(:,:,model.T == t)));
        
        
        % Set axes options
        ax(I).YDir = 'normal';
        ax(I).CLim = [0 9];
        ax(I).FontSize = 14;

        ax(I).LineWidth = 2;
        
        if j > 1
            ax(I).YTick = [];
        else
            ax(I).YTick = 1:4:ax(I).YLim(2);
            ax(I).YTickLabel = cellfun(@(x)sprintf('10^{%.2g}',log10(model.CFU(x))),num2cell(ax(I).YTick),'UniformOutput',false);
        end
        if i < n
            ax(I).XTick = [];
        else
            ax(I).XTick = 1:4:ax(I).XLim(2);
            ax(I).XTickLabel = cellfun(@(x)sprintf('10^{%.2g}',log10(model.PFU(x))),num2cell(ax(I).XTick),'UniformOutput',false);
%             ax(I).XTick(1)   = ax(I).XTick(1)+0.5;
%             ax(I).XTick(end) = ax(I).XTick(end)-0.5;
        end
        
        if i == 1
            title(sprintf('l = %d',1e4/nGrid(j)),'FontWeight','Normal')
        end
    end
end

% Add labels to graph
ax0 = axes('Position',[0 0 1 1],'Visible','off');

% X Label
hx = text(0.5,0.025,'Initial phage density (PFU / ml)');
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
hcbar = text(0.97,0.515,'Bacterial density at time T = 15 h (CFU / ml)');
hcbar.Rotation = 90;
hcbar.FontSize = 14;
% hy.FontWeight = 'bold';
hcbar.HorizontalAlignment = 'center';

% A/B/C/D Label
hA = text(mar+1.2*sep,3*(hei+sep)+hei+mar+0.5*sep,'d','Units','Pixels');
hA.Color = 'k';
hA.FontWeight = 'bold';
hA.FontSize = 16;

% hA = text(mar+1.2*sep,3*(hei+sep)-sep+hei+mar+0.8*sep,'{\Delta}T = 1x','Units','Pixels');
% hA.Color = 'k';
% hA.FontWeight = 'bold';
% hA.FontSize = 10;

hB = text(mar+1.2*sep,2*(hei+sep)+hei+mar+0.5*sep,'d','Units','Pixels');
hB.Color = 'k';
hB.FontWeight = 'bold';
hB.FontSize = 16;

hB = text(mar+1.2*sep,2*(hei+sep)-sep+hei+mar+0.9*sep,'{\Delta}T = 0.5x','Units','Pixels');
hB.Color = 'k';
hB.FontWeight = 'bold';
hB.FontSize = 10;

hC = text(mar+1.2*sep,(hei+sep)+hei+mar+0.5*sep,'d','Units','Pixels');
hC.Color = 'k';
hC.FontWeight = 'bold';
hC.FontSize = 16;

hC = text(mar+1.2*sep,(hei+sep)-sep+hei+mar+0.9*sep,'{\Delta}T = 0.25x','Units','Pixels');
hC.Color = 'k';
hC.FontWeight = 'bold';
hC.FontSize = 10;

hD = text(mar+1.2*sep,hei+mar+0.5*sep,'d','Units','Pixels');
hD.Color = 'k';
hD.FontWeight = 'bold';
hD.FontSize = 16;

hD = text(mar+1.2*sep,-sep+hei+mar+0.9*sep,'{\Delta}T = 0.1x','Units','Pixels');
hD.Color = 'k';
hD.FontWeight = 'bold';
hD.FontSize = 10;

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


if ~exist('../Fig_S1','dir')
    mkdir('../Fig_S1')
end
print('../Fig_S1/FigS1.tif','-dtiff')