clearvars; close all;

path = strrep(pwd,'analysis','cpp/data');

% Define functions

% Sheilding function 1
S1 = @(C,i) max(0, 1-C*i/C^(2/3));

% Sheilding function 2
S2 = @(C,i) max(0, double(C<64) * (1-C*i/C) + double(C>64) * (1-C*i/(4*C^(2/3))));

% Sheilding function 3
S3 = @(C,i,k) max(0,min(1-i, 1 - C*i/(C - (C^(1/3)-k)^3)));

% Sheilding function 4
S4 = @(C,i,zeta) max(0,min(1-i,exp(-zeta*(1.3)^(1/3)*( C^(1/3) - (C-C*i)^(1/3) ))));

% Diffusion limited
ColonySize = [10 32 100 316 1000];
InfectionFraction = [0.2 0.4 0.6 0.8];
gamma = [10^2 10^2.5 10^3 10^3.5 10^4 10^4.5 10^5 10^5.5 10^6];

% Make mesh grid of sample points
[i,c] = meshgrid(InfectionFraction,ColonySize);

% Load he data
for g = 1:9

    switch g
        case 1
            filePath = [path '/PhageInfectionProbability/gamma_1e2.txt'];
        case 2
            filePath = [path '/PhageInfectionProbability/gamma_1e2_5.txt'];
        case 3
            filePath = [path '/PhageInfectionProbability/gamma_1e3.txt'];
        case 4
            filePath = [path '/PhageInfectionProbability/gamma_1e3_5.txt'];
        case 5
            filePath = [path '/PhageInfectionProbability/gamma_1e4.txt'];
        case 6
            filePath = [path '/PhageInfectionProbability/gamma_1e4_5.txt'];
        case 7
            filePath = [path '/PhageInfectionProbability/gamma_1e5.txt'];
        case 8
            filePath = [path '/PhageInfectionProbability/gamma_1e5_5.txt'];
        case 9
            filePath = [path '/PhageInfectionProbability/gamma_1e6.txt'];
    end

    SuccessRate{g} = nan(5,4);

    fh = fopen(filePath);
    for l = 1:5
        fgetl(fh);
        fgetl(fh);
        for m = 1:4
            fgetl(fh);
            fgetl(fh);
            SuccessRate{g}(l,m) = str2double(strrep(fgetl(fh),'	The probability of hitting a non-infected cell is: ',''));
            fgetl(fh);
        end
        fgetl(fh);
        fgetl(fh);
    end
    fclose(fh);

    % Model 3 and 4 has fit parameter, so we define chi2
    chi2_k{g}    = @(k)   sum(sum((SuccessRate{g}-arrayfun(@(c,i)S3(c,i,k),c,i)).^2));
    chi2_zeta{g} = @(zeta)sum(sum((SuccessRate{g}-arrayfun(@(c,i)S4(c,i,zeta),c,i)).^2));


end

load Adsorption.mat %% This data comes from a previous published paper (see repository: https://github.com/RasmusSkytte/GrowingMicrocolonyProtection for details)

k = arrayfun(@(n)fminsearch(chi2_k{n},1),1:9);
zeta = arrayfun(@(n)fminsearch(chi2_zeta{n},1),1:9);


fh = figure(1); clf;
fh.Position = [1725 -79 1113 488];
s = semilogx(v_eta(end:-1:1),1./zeta,'o-','LineWidth',2);
s.MarkerFaceColor = s.Color;
xlabel('\eta ({\mu}m^3/h)')
ylabel('\zeta')
hold off
xlim([10^(2.5) 10^(5.5)])
ylim([0.5 3.5])
grid on;
ax = gca;
ax.YTick = 0.5:0.5:3.5;
ax.FontSize = 20;

ax.LineWidth = 2;

% legend show;
drawnow;

if ~exist('Fig_S3','dir')
    mkdir('Fig_S3')
end
print('Fig_S3/zeta_vs_eta.tif','-dtiff')


% Run over data
for g = 5

    % For each gamma, fit the functions to the data
    fh = figure(2); clf;
    fh.Position = [1500  -150   750   600];

    cc = logspace(0,log10(ColonySize(end)),100);
    ii = linspace(0,InfectionFraction(end),100);
    [CC,II] = meshgrid(cc,ii);

    fa = 0;%0.65;

    ax = axes;
    ax.NextPlot = 'add';
    view([-0.5 1 0.5])
    %     shading flat;
    set(ax,'XScale','log')
    ax.FontSize = 30;
    ax.Box = 'on';
    ax.XTick = [1 10 100 1000];
    ax.YTick = [0 0.4 0.8];
    ax.ZTick = [0 0.5 1.0];
    ax.YTickLabel = {'0.0', '0.4', '0.8'};
    ax.ZTickLabel = {'0.0', '0.5', '1.0'};
    xlim([1 1000])
    ylim([0 0.8])
    zlim([0 1])
    xl = xlabel('B + I');
    xl.Rotation = -12;
    yl = ylabel('I / B');
    yl.Rotation = 48;
    zlabel('S')

    fh.Color = [1 1 1];

    ax.LineWidth = 2;
    
    for n = 1:4
        for p = 1:numel(SuccessRate{g})

            switch n
                case 1
                    m = S1(c(p),i(p));
                case 2
                    m = S2(c(p),i(p));
                case 3
                    m = S3(c(p),i(p),k(g));
                case 4
                    m = S4(c(p),i(p),zeta(g));
            end

            if m > SuccessRate{g}(p)
                s(p) = scatter3(ax,c(p),i(p),SuccessRate{g}(p),50,'r^','filled');
            else
                s(p) = scatter3(ax,c(p),i(p),SuccessRate{g}(p),50,'rv','filled');
            end

        end

        switch n
            case 1
                s(p+1) = mesh(ax,cc,ii,arrayfun(S1,CC,II),'FaceAlpha',fa);
                string = 'A';
            case 2
                s(p+1) = mesh(ax,cc,ii,arrayfun(S2,CC,II),'FaceAlpha',fa);
                string = 'B';
            case 3
                s(p+1) = mesh(ax,cc,ii,arrayfun(@(c,i)S3(c,i,k(g)),CC,II),'FaceAlpha',fa);
                string = 'C';
            case 4
                s(p+1) = mesh(ax,cc,ii,arrayfun(@(c,i)S4(c,i,zeta(g)),CC,II),'FaceAlpha',fa);
                string = 'D';
        end
        ax0 = axes();
        ax0.Visible = 'off';
        text(ax0,-0.1, 1,string,'FontSize',20,'FontWeight','bold','HorizontalAlignment','center')


        shading(ax,'flat')

        set(gcf,'InvertHardcopy','off')

        if ~exist('Fig_S2','dir')
            mkdir('Fig_S2')
        end
        print(sprintf('Fig_S2/ShieldingFunction_S%d_eta_%.1g.tif',n,gamma(g)),'-dtiff')

        delete(s)
        delete(ax0)
    end



end
