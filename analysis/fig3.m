clearvars; close all;

path = strrep(pwd,'analysis','cpp/data');

% Choose which point to show dynamics for:
B0 = 1e4;
P0 = 1e5;

% Create figures
fh = figure();
fh.Position = [1600 -100 750 400];

labels = {'A','B','C','D'};

for m = [1 2 3 4]
    % Reset
    clf;
    
    % Load data
    subpaths = ls(sprintf('%s/Model_%d',path,m));
    subpaths = strsplit(subpaths(1:end-1),' ');
    importPath = sprintf('%s/Model_%d/%s/data_Model_%d.mat',path,m,subpaths{end},m);
    model = importdata(importPath);
    
    % Prepare axes
    ax = gca;
    ax.Box = 'on';
    ax.XLabel.String = 'Time (hours)';
    ax.YLabel.String = 'Population Size';
    ax.FontSize = 28;
    ax.YScale = 'log';
    
    ax.XLim = [0 15];
    ax.YLim = [1 1e12];
    
    ax.YTick = 10.^(0:4:12);
    
    ax.LineWidth = 2;
    
    % Add label
    h = text(0.5,10^11,labels{m});
    h.FontSize = 28;
    h.FontWeight = 'bold';

    grid on;
    ax.GridAlpha = 0.5;
    
    % Plot data
    hold on;
        
    plot(model.T,max(1e-3,squeeze(model.B(model.CFU==B0, model.PFU==P0,:))),'b-','LineWidth',3)
    plot(model.T,max(1e-3,squeeze(model.P(model.CFU==B0, model.PFU==P0,:))),'r:','LineWidth',3)
    
    hold off;
    drawnow;
    if ~exist('Fig_3','dir')
        mkdir('Fig_3')
    end
    print('-r300',sprintf('Fig_3/Dynamics_Model_%d.tif',m),'-dtiff')

    
end
close(fh);