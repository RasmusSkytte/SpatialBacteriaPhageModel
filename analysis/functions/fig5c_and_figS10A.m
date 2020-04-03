clearvars; close all;

path = strrep(pwd,'analysis/functions','cpp/data');

NUM = 75;

threshold = 4*pi/3*50^3;

path = [path '/Experiment_Conditions/beta_400_delta_0.003_eta_13200_tau_1.00_D_P_3000_nGrid_100_boundaryModifier_10'];

% Define the points to sample over:
T_i = 0:8;

% Load experimental survival curve:
CFU = ...
    [ 81 0  0  0;
    66 0  0  0;
    76 0  0  0;
    68 0  0  0;
    84 13 10 4;
    82 31 15 25;
    66 64 56 51;
    61 76 54 NaN;
    96 65 84 78];

t_d = (0:8)';

m_c = nanmean(CFU(:,1));
s_c = nanstd(CFU(:,1))/sqrt(sum(~isnan(CFU(:,1))));

m_p = nanmean(CFU(:,2:4),2);
s_p = nanstd(CFU(:,2:4),[],2)./sqrt(sum(~isnan(CFU(:,2:4)),2));

f_d = m_p/m_c;
s_d = sqrt((s_p./m_p).^2 + (s_c/m_c)^2).*f_d;
s_d(isnan(s_d)) = 0;

fh = figure(2); clf;
fh.Position = [ 1520 159 1437 600];

ax1 = axes(); box on;

ax1.XLim = [0 8];
ax1.XTick = [];
ax1.YLim = [0 1.1];
ax1.YTick = 0.0:0.5:1;
ax1.YTickLabel{1} = '0.0';
ax1.YTickLabel{3} = '1.0';

ax1.YLabel.String = 'Survival Fraction';
tx = text(ax1,-0.5, 1.05,'c','FontSize',20,'FontWeight','bold','HorizontalAlignment','center');

ax1.Position = [0.07    0.12    0.92    0.875];
ax1.FontSize = 20;

ax1.XLabel.String = 'Time of phage addition (h)';
ax1.XLim = [0 8];
ax1.XTick = 0:2:8;

ax1.NextPlot = 'add';

ax1.LineWidth = 2;


h1 = errorbar(ax1,t_d,f_d,s_d,'bo','MarkerSize',8,'CapSize',8,'MarkerFaceColor','Auto','LineWidth',2); hold on;
% Allocate array
nutrient = nan(length(T_i),1e3,1e4);    % Nutrient consumed
biomass  = nan(length(T_i),1e3,1e4);    % Biomass

% Loop over zeta
z = [-1 0.4 0.2 0.1];
for i = 1:length(z)
    
    % Loop over times to store the data
    for t = 1:length(T_i)
        
        if z(i) > 0
            testPath = [path '_full_zeta_' num2str(z(i))];
        elseif z(i) == 0
            testPath = [path '_full_noSheilding'];
        elseif z(i) == -1
            testPath = [path '_Model_3'];
        end
        
        if exist(testPath,'dir')
            
            % Construct full path
            fullPath = [testPath sprintf('/T_i_%.1f/',T_i(t))];
            
            % Loop over elements in the folder
            d = dir(fullPath);
            d(strcmpi({d.name},'.')) = [];
            d(strcmpi({d.name},'..')) = [];
            %                 d(strcmpi({d.name},'-1')) = [];
            [~,I] = sort(cellfun(@(d)str2double(d),{d.name}));
            d = d(I);
            
            for n = 1:length(d)
                if strcmpi(d(n).name,'.') || strcmpi(d(n).name,'..')
                    continue
                else
                    k = str2double(d(n).name);
                    
                    if k == -1
                        continue
                    end
                    if k >= NUM
                        continue
                    end
                    
                    if ~exist(sprintf('%s/%d/Completed.txt',fullPath,k),'file')
                        warning('Not Ready: zeta = %.3f, T_i = %.f h: run # %d\n',z(i),T_i(t),k)
                        continue;
                    end
                    
                    %                             fprintf('zeta = %.3f, T_i = %.f h: run # %d\n',z,T_i(t),k)
                    
                    lh = fopen(sprintf('%s/%d/log.txt',fullPath,k));
                    try
                        while ~feof(lh)
                            str = [fgetl(lh) ';'];
                            eval(str);
                        end
                    catch
                    end
                    fclose(lh);
                    
                    kk = k;
                    
                    pop = importdata(sprintf('%s/%d/PopulationSize.txt',fullPath,k));
                    
                    if k >= 0
                        if size(pop.data,2) > 5
                            nutrient(t,kk+1,1:size(pop.data,1)) = pop.data(:,6);
                            nutrient(t,kk+1,(size(pop.data,1)+1):(nSamp*T_end+1)) = pop.data(end,6); % Repeat last
                            
                            biomass(t,kk+1,1:size(pop.data,1)) = sum(pop.data(:,2:3),2);
                            biomass(t,kk+1,(size(pop.data,1)+1):(nSamp*T_end+1)) = sum(pop.data(end,2:3),2); % Repeat last
                        end
                    end
                end
            end
        end
    end
    
    m_r = nan(length(T_i),1);
    e_r = nan(length(T_i),1);
    
    % Take average of each experiment
    for t = 1:length(T_i)
        tt = (T_i(t)+16)*nSamp+1; % Extract at time T_i + 10
        
        % Check if run is completed
        I = ~isnan(biomass(t,:,tt));
        
        m_r(t) = mean(nutrient(t,I,tt)>threshold);
        e_r(t) = std(nutrient(t,I,tt)>threshold);
        
        e_r(t) = e_r(t)/sqrt(sum(I));
        
    end
    
    % Define plotting style
    cc = lines(length(z));
    color = cc(i,:);
    angle = 90/(length(z)-1)*(i-1)+90;
    switch i
        case 1
            style = ':';
            color = 'k';
        case 2
            style = '--';
            color = 'r';
        case 3
            style = '-.';
            color = [0 0.6 0];
        case 4
            style = '-';
        case 5
            style = '.';
    end
    
    % Compare with experimental data
    s = shadedErrorBar(T_i', m_r, e_r,'lineProps',{style,'Color',color,'LineWidth',4});
    s.edge(1).LineWidth = 2;
    s.edge(2).LineWidth = 2;
    
    h = hatchfill(s.patch,'single',angle,10);
    h.Color = color;
    h.LineWidth = 2;
    delete(s.patch);
    
    
    uistack(h1,'top');
    
    drawnow;
    
end

drawnow;


if ~exist('../Fig_5','dir')
    mkdir('../Fig_5')
end

print('-r400','../Fig_5/lower.tif','-dtiff')

% Change the label
tx.String = 'D_P = 3000 {\mu}m^2/h';
tx.Position(1) = 1.0;
tx.Position(2) = 1.0;


if ~exist('../Fig_S10','dir')
    mkdir('../Fig_S10')
end
print('-r400','../Fig_S10/upper.tif','-dtiff')

