path = strrep(pwd,'analysis','cpp/data');
extractData(path)

function extractData(varargin)

if nargin == 0  % No path as input: Set current path as master
    
    % Get current directory
    d = dir;
    path = '';
    
else
    
    % Check varargin has valid path
    if exist(varargin{1},'dir')
        d = dir(varargin{1});
    else
        return;
    end
    
    % Set varargin as master path
    path = [varargin{1} '/'];
    
end

% The path does not contain any children
if length(d) < 3
    return
end

% Check if the folder contains the right format
str = '';
for i = 1:length(d)
    if strcmpi(d(i).name,'.') || strcmpi(d(i).name,'..')
        continue
    end
    
    % We do not want images
    if ~isempty(strfind(d(i).name,'.fig'))
        continue;
    end
    if ~isempty(strfind(d(i).name,'.png'))
        continue;
    end
    
    % If it is a folder
    if d(i).isdir
        str = regexp(d(i).name,'_','split');
        
        % Does the folder have the correct name?
        if strcmpi(str{1},'CFU') && strcmpi(str{3},'PFU')
            
            
            % Get the values searched over
            CFU = [];
            PFU = [];
            
            for i = 3:length(d)
                
                if ~isempty(strfind(d(i).name,'.fig'))
                    continue;
                end
                if ~isempty(strfind(d(i).name,'.png'))
                    continue;
                end
                if ~isempty(strfind(d(i).name,'.mat'))
                    continue;
                end
                
                str = regexp(d(i).name,'_','split');
                if numel(str) ~= 4
                    continue
                end
                CFU = [CFU eval(strrep(str{2},'1e','10^'))];
                PFU = [PFU eval(strrep(str{4},'1e','10^'))];
            end
            
            CFU = sort(unique(CFU));
            PFU = sort(unique(PFU));
                       
            if length(CFU) == 1 || length(PFU) == 1
                return
            end
            
            % Check if data is already analyzed
            if exist([path 'data.mat'],'file')
                cfu = CFU;
                pfu = PFU;
                load([path 'data.mat']);
                if length(cfu) == length(CFU) && all(cfu == CFU) && length(pfu) == length(PFU) && all(pfu == PFU)
                    return;
                end
                CFU = cfu;
                PFU = pfu;
            end
            
            f = nan(length(CFU),length(PFU),5000);
            n = nan(length(CFU),length(PFU),5000);
            b = nan(length(CFU),length(PFU),5000);
            B = nan(length(CFU),length(PFU),5000);
            C = nan(length(CFU),length(PFU),5000);
            I = nan(length(CFU),length(PFU),5000);
            D = nan(length(CFU),length(PFU),5000);
            P = nan(length(CFU),length(PFU),5000);
            nC = nan(length(CFU),length(PFU),5000);
            foo = nan(length(CFU),length(PFU),5000);
            
            for c = 1:length(CFU)
                for p = 1:length(PFU)
                    
                    if ~exist([path sprintf('CFU_1e%.2f_PFU_1e%.2f/Completed.txt',log10(CFU(c)),log10(PFU(p)))],'file')
                        return
                    end
                    
                    display([path sprintf('CFU_1e%.2f_PFU_1e%.2f/',log10(CFU(c)),log10(PFU(p)))])
                    
                    fh = fopen([path sprintf('CFU_1e%.2f_PFU_1e%.2f/log.txt',log10(CFU(c)),log10(PFU(p)))]);
                    try
                        while ~feof(fh)
                            eval([fgetl(fh) ';']);
                        end
                    catch
                    end
                    fclose(fh);
                    
                    pop = importdata([path sprintf('CFU_1e%.2f_PFU_1e%.2f/PopulationSize.txt',log10(CFU(c)),log10(PFU(p)))]);
                    pop.data = pop.data(pop.data(:,1) <= 20,:);
                    T_end = 20;
                       
                    % If the fast exit is used, repeat the last data point
                    if size(pop.data,1) < (T_end*nSamp+1)
                        % Rows remaining
                        k = (T_end*nSamp+1) - size(pop.data,1);

                        % Repeat the last data point
                        pop.data = [pop.data; repmat(pop.data(end,:),k,1)];
                    end

                    B(c,p,1:(T_end*nSamp+1)) = pop.data(:,2);
%                     I(c,p,1:(T_end*nSamp+1)) = pop.data(:,3);
                    P(c,p,1:(T_end*nSamp+1)) = pop.data(:,4);
                    f(c,p,1:(T_end*nSamp+1)) = pop.data(:,5);
                    n(c,p,1:(T_end*nSamp+1)) = pop.data(:,6);
                    b(c,p,1:(T_end*nSamp+1)) = sum(pop.data(:,2:3),2);

                end
            end
            

            % Save the data
            T = (0:(T_end*nSamp))/nSamp;
            B = B(:,:,1:(T_end*nSamp+1));
%             I = I(:,:,1:(T_end*nSamp+1));
            P = P(:,:,1:(T_end*nSamp+1));
            f = f(:,:,1:(T_end*nSamp+1));
            n = n(:,:,1:(T_end*nSamp+1));
            b = b(:,:,1:(T_end*nSamp+1));
            
            foo = strsplit(path,'/');
            save([path 'data_' foo{end-2} '.mat'],'PFU','CFU','T','n','b','f','B','I','P');
            
            return % Do not go deeper upon success.
        else % Continue into folder
            extractData([path d(i).name])
        end
    end
end


if isempty(str)
    return;
end

end