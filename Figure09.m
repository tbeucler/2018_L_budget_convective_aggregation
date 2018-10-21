%% Figure09.m
% tbeucler - 9/23/2018
% tbeucler - 8/20/2018
% Plot length scale rates versus sea surface temperature
% for both radiation schemes

close all; fclose('all'); clearvars;

%% 1. Parameters and initialization of figure and table
% User's choice
Lmet = 'LAM'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Physical parameters
spd = 24*3600; % Number of seconds per day

% Model parameters
i2 = 25; % Beginning of second day of simulation
i31 = 24*30+1; % Beginning of thirty-first day of simulation
i75 = 24*75; % End of the seventy-fifth day of simulation
gro = i2:(i31-1); mat = i31:i75; % First month and rest of the simulation
SST_array = 280:5:310; NSST = numel(SST_array); % Sea surface temperatures [K]
rad_array = {'cam','rrtm'}; Nrad = numel(rad_array); % Radiation scheme
toplot = struct; % Structure for reduced growth and mature rates

% Figure parameters
ANOT = {'(a) CAM','(b) RRTM'}; % Panel annotations
colors_rgb = [117 112 179; 217 150 0; 27 158 119; 231 41 138]; % LW/SW/SEF/ADV
cmap = colors_rgb./255; % Colorbar from Tim & Allison
DelX = 0.1; DelY = 0.05; % Legend dimensions
f_fig = {'lw','sw','sef','adv'}; Nf_fig = numel(f_fig); % Fields plotted in figure
fz = 12; % Fontsize
lw = 2.5; % Linewidth
sz = 75; % Marker size
LEG = {'LWgrowth','SWgrowth','SEFgrowth','ADVgrowth','LWmature',...
    'SWmature','SEFmature','ADVmature'}; % Legend for L rates
LOC='Northoutside'; OR='horizontal'; % Location and orientation of colorbar
XLAB = 'Sea Surface Temperature [K]'; % x label
YLAB = '$\dot{L_{i}}/L\ \left[1/\mathrm{day}\right]$'; % y label
if strcmp(Lmet,'Km1')==1, YLIM = [-0.6 0.75]; % y bounds of figure
elseif strcmp(Lmet,'LAM')==1, YLIM = [-0.35 0.35]; % y bounds of figure
end

for irad = 1:Nrad, rad = rad_array{irad};
    for iSST = 1:NSST, SST = SST_array(iSST);
        
        % Name and load each experiment
        if irad==3
            if iSST==1, name = 'mdradsfc'; % MD RAD+SFC experiment
            elseif iSST==2, name = 'mdsfc'; % MD RAD experiment
            elseif iSST==3, name = 'mdrad'; % MD SFC experiment
            elseif iSST==4, name = 'sqcam'; % CAM square experiment
            elseif iSST==5, name = 'sqrrtm'; % RRTM square experiment
            elseif iSST==6, name = 'bsqrrtm'; % Big square experiment
            end
        else, name = [rad,num2str(SST)];
        end; load(['MAT_DATA',filesep,name,'.mat']); % Load data
        
        %% 2. Prepare reduced data for figure
        
        f = fieldnames(dat.(Lmet)); Nf = numel(f);
        for i = 1:Nf_fig
            % Length scale tendencies averaged over Gro and Mat phases
            toplot.Gro.(f_fig{i})(iSST,irad) = spd*squeeze(mean(...
                dat.(Lmet).(f_fig{i})(:,:,gro)./dat.(Lmet).mse(:,:,gro),3));
            toplot.Mat.(f_fig{i})(iSST,irad) = spd*squeeze(mean(...
                dat.(Lmet).(f_fig{i})(:,:,mat)./dat.(Lmet).mse(:,:,mat),3));
            % 3 x Standard mean error of length scale tendency
            % (STD/sqrt(Sample size))
            toplot.MEg.(f_fig{i})(iSST,irad) = 3*spd*squeeze(std(dat.(Lmet).(...
                f_fig{i})(:,:,gro)./dat.(Lmet).mse(:,:,gro),0,3))/sqrt(numel(gro));
            toplot.MEm.(f_fig{i})(iSST,irad) = 3*spd*squeeze(std(...
                dat.(Lmet).(f_fig{i})(:,:,mat)./dat.(Lmet).mse(:,:,mat),0,3))/...
                sqrt(size(dat.(Lmet).(f_fig{i})(:,:,mat),3));
        end
    end
end

%% 3. Figure

figure; set(gcf,'Position',[10 10 1200 600]); % Figure dimensions
for irad = 1:Nrad, S(irad) = subplot(Nrad,1,irad); box on; % Subplot
    line([SST_array(1)-2.5 SST_array(end)+2.5],[0 0],'Linewidth',lw/2,...
        'color',[0.6 0.6 0.6]); hold on; % 0 line for reference
    for i = 1:Nf_fig
        for iSST = 1:NSST, SST = SST_array(iSST);
            errorbar(SST_array+(i-2.5)/2,toplot.Gro.(f_fig{i})(:,irad),...
                toplot.MEg.(f_fig{i})(:,irad),'Linewidth',lw,'Linestyle','none',...
                'color',cmap(i,:)); hold on; % Error bar for mean square error
            errorbar(SST_array+(i-2.5)/2,toplot.Mat.(f_fig{i})(:,irad),...
                toplot.MEm.(f_fig{i})(:,irad),'Linewidth',lw,'Linestyle','none',...
                'color',cmap(i,:)); hold on; % Error bar for mean square error
            G(i) = scatter(SST_array+(i-2.5)/2,toplot.Gro.(f_fig{i})(:,irad),...
                'MarkerFaceColor','w','MarkerEdgeColor',cmap(i,:),...
                'SizeData',sz,'Linewidth',lw); hold on; % Average L rate for month 1
            M(i) = scatter(SST_array+(i-2.5)/2,toplot.Mat.(f_fig{i})(:,irad),...
                'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),...
                'SizeData',sz); hold on; % Average L rate for later months
        end
    end; grid on; xlim([SST_array(1)-2.5 SST_array(end)+2.5]); ylim(YLIM); % x and y limits
    ylabel(YLAB,'Interpreter','Latex','Fontsize',fz); % y label
    set(gca,'TickLabelInterpreter','Latex','Fontsize',fz); % Axes appearance
end; xlabel(XLAB,'Interpreter','Latex','Fontsize',fz); % x label

%% 4. Legend

% Adjust legend position
LOBJ = legend([G M],LEG); G = gca; POS = G.Position;
set(LOBJ,'Location',LOC,'orientation',OR,'Fontsize',fz,'Interpreter','Latex');
set(gca,'Position',POS); LOBJ_X = LOBJ.Position(3);
LOBJ.Position(1) = (1-LOBJ_X)/2; LOBJ.Position(3) = LOBJ_X;

% Optimize space between subplots
set(S(1),'XTickLabel',''); LOBJ.Position(2) = 0.95;
S(2).Position = [0.1 0.1 0.8 0.4]; S(1).Position = [0.1 0.525 0.8 0.4];

% (a)/(b) Annotations for each subplot
for irad = 1:2
    A(irad) = annotation('textbox',...
        [S(irad).Position(1)+S(irad).Position(3)-DelX ...
        S(irad).Position(2)+S(irad).Position(4)-DelY ...
        DelX DelY],'Fontsize',fz+3,'String',ANOT{irad},...
        'BackgroundColor',[1 1 1],'Interpreter','Latex','Linewidth',lw,...
        'HorizontalAlignment','center'); hold on;
end

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'Figure09';
elseif strcmp(Lmet,'Km1')==1, name = 'FigureS9';
end; gcfsavepdf([basedir,'PDF_DATA',filesep,name,'.pdf']); % Save plot
