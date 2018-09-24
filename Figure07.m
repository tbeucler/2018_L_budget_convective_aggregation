%% Figure07.m
% tbeucler - 9/23/2018
% tbeucler - 8/30/2018
% Decomposes the radiative length scale tendency by plotting
% time-averaged longwave and shortwave tendencies for month 1 & after
% in (L,AGG) space: Arrows depict stretching/shrinking from L to L_DIAB

close all; fclose('all'); clearvars;

%% 1. Parameters
% User's choice
Lmet = 'LAM'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Physical parameters
spd = 24*3600; % Number of seconds per day

% Data parameters
i2 = 25; % Beginning of second day of simulation
i31 = 24*30+1; % Beginning of thirty-first day of simulation
i75 = 24*75; % End of the seventy-fifth day of simulation
gro = i2:(i31-1); mat = i31:i75; % First month and rest of the simulation
cas_array = {'cam300','cam280','cam310','rrtm300','mdsfc'}; Ncas = numel(cas_array);

% Figure parameters
cmap = parula(7); % Parula color scale if parallel with figure 6 of Wing & Cronin
cmap(7,:) = [139 69 19]/255; % Make the yellow saddlebrown
COL = cat(1,cmap(5,:),cmap(1,:),cmap(7,:),cmap(5,:),[1 0 0]); % Color scale
f_rad = {'lw','sw'}; Nf_rad = numel(f_rad); % Radiation fields
fz = 12; % Fontsize
LEGarray = {'$\mathrm{CAM\ }300\mathrm{K}$','$\mathrm{CAM\ }280\mathrm{K}$',...
    '$\mathrm{CAM\ }310\mathrm{K}$',...
    '$\mathrm{RRTM\ }300\mathrm{K}$','$\mathrm{MD\ SFC}$'}; % Case legends
lw = 2; % Linewidth
sz = 50; % Marker size
TITarray = {'$\mathrm{\left(a\right)\ Longwave-Growth\ phase}$',...
    '$\mathrm{\left(b\right)\ Longwave-Mature\ phase}$',...
    '$\mathrm{\left(c\right)\ Shortwave-Growth\ phase}$',...
    '$\mathrm{\left(d\right)\ Shortwave-Mature\ phase}$'}; % Title array
XLAB = '${\cal L}_{i}\ \left[\mathrm{1000km}\right]$'; % X label
YLAB = ['$\left\langle \dot{\varphi}_{i}\right\rangle /\left\langle \',...
    'varphi\right\rangle \ \left[1/\mathrm{day}\right]$']; % Y label
if strcmp(Lmet,'Km1')==1
Xmin = [[0 0];[0 0]]; Xmax = [[2.25 2.25];[4.5 4.5]]; % X bounds
Ymin = [[0 0];[0 0]]; Ymax = [[0.35 0.66];[0.35 0.66]]; % Y bounds
elseif strcmp(Lmet,'LAM')==1
Xmin = [[0 0];[0 0]]; Xmax = [[4.5 4.5];[6.1 6.1]]; % X bounds
Ymin = [[0 0];[0 0]]; Ymax = [[0.35 0.64];[0.35 0.64]]; % Y bounds
end

%% 2. Prepare reduced data for figure
figure; set(gcf,'Position',[50 50 750 750]); % Figure's position

for i = 1:Ncas, load(['MAT_DATA',filesep,cas_array{i},'.mat']); % Load data
    
    % Average data in growth and mature phases
    for ipha = 1:2 % First month [1] Last months [2] of simulation
        if ipha==1, dom = i2:i31; else, dom = i31:i75; end
        for irad = 1:Nf_rad, L = nanmean(dat.(Lmet).mse(:,:,dom),3); % Here bc units
            
            % Decomposition of length scale tendency
            AGG = nanmean(dat.AGG.(f_rad{irad})(:,:,dom),3); % Aggregation rate [1/s]
            if strcmp(Lmet,'LAM')==1
                LLW = nanmean(dat.(Lmet).(f_rad{irad})(:,:,dom),3); % L tendency [m/s]
                L_DIAB = L+LLW/AGG; % Diabatic length scale from agg. and L. rates [m]
            elseif strcmp(Lmet,'Km1')==1
                LLW = nanmean(dat.(Lmet).(f_rad{irad})(:,:,dom)./...
                    (dat.(Lmet).mse(:,:,dom)),3); % L rate [1/s]
                L_DIAB = L/(1-LLW/AGG); % Diabatic length scale [m]
                L_DIAB(L_DIAB>1e6*Xmax(ipha,irad)) = 1e6*Xmax(ipha,irad); % Graph cannot fit diab. L larger than Xmax
                L_DIAB(L_DIAB<0) = 1e6*Xmax(ipha,irad); % Set super shrinking diab. L to Xmax
            end
            
            % Convert rates and scales to [1/day] and [1000 km] for plotting
            AGG = spd*AGG; L_DIAB = L_DIAB/1e6; L = L/1e6;
            
            %% 3. Figure
            S(ipha,irad) = subplot(2,2,ipha+2*(irad-1)); % Subplot
            
            % Plot arrows and lines
            FROM(i,ipha,irad) = scatter(L,AGG,'MarkerEdgeColor',COL(i,:),...
                'SizeData',sz,'Linewidth',lw); hold on; % Arrow from FROM
            LINE(i,ipha,irad) = line([min(L_DIAB,L) max(L_DIAB,L)],[AGG AGG],...
                'Linewidth',lw,'Color',COL(i,:)); hold on; % Arrow's line
            TO(i,ipha,irad) = scatter(L_DIAB,AGG,'MarkerEdgeColor',COL(i,:),...
                'SizeData',sz,'Linewidth',lw); hold on; % Arrow to TO
            if i~=4, set([FROM(i,ipha,irad) TO(i,ipha,irad)],'MarkerFaceColor',COL(i,:));
            else, set(LINE(i,ipha,irad),'Linestyle',':'); % Dot RRTM arrow
            end; grid on;
            if L_DIAB>L, set(TO(i,ipha,irad),'Marker','>'); % Stretching arrow
            else, set(TO(i,ipha,irad),'Marker','<'); % Shrinking arrow
            end
            
            % Axes
            title(TITarray{ipha+2*(irad-1)},'Fontsize',fz,'Interpreter','Latex');
            if ipha==1, ylabel(YLAB,'Fontsize',fz,'Interpreter','Latex'); end
            if irad==2, xlabel(XLAB,'Fontsize',fz,'Interpreter','Latex'); end
            xlim([Xmin(ipha,irad) Xmax(ipha,irad)]);
            ylim([Ymin(ipha,irad) Ymax(ipha,irad)]);
            
        end
    end
end

%% 4. Legend
for ipha = 1:2
    for irad = 1:Nf_rad
        
        set(S(ipha,irad),'Fontsize',fz,'Linewidth',lw,'TickLabelInterpreter','Latex');
        if ipha+2*(irad-1)==1, LEG = legend(TO(:,ipha,irad),LEGarray,...
                'Fontsize',fz,'Interpreter','Latex','orientation','horizontal');
        end % Write legend and adjust subplot position
        S(ipha,irad).Position(2) = S(ipha,irad).Position(2)-0.05*(3-irad);
        
    end
end
LEG.Position(1) = 0.5*(-LEG.Position(3)+S(1,1).Position(1)+...
    S(2,1).Position(1)+S(2,1).Position(3)); % Set horizontal legend position
LEG.Position(2) = S(1,1).Position(2)+1.2*S(1,1).Position(4); % Set vert. leg. position

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'Figure07';
elseif strcmp(Lmet,'Km1')==1, name = 'FigureS7';
end; gcfsavepdf([basedir,'PDF_DATA',filesep,name,'.pdf']); % Save plot
