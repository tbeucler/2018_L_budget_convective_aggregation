%% Figure05.m
% tbeucler - 9/23/2018
% tbeucler - 8/30/2018
% For the reference simulation LC300CAM plots
% (a) The length scale tendencies [m/s] (b) the aggregation rates [1/d]
% versus simulation time [day]

fclose('all'); close all; clearvars;

%% 1. Parameters and data
% User's choice
Lmet = 'LAM'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Physical parameters
DAY = 24; % Number of timesteps in 1 day
WEEK = 7*DAY; % Number of timesteps in 1 week
spd = 24*3600; % Number of seconds per day

% Data parameters
i2 = 25; % First timestep of the second day
i31 = 24*30+1; % Beginning of thirty-first day of simulation
SST = 300; % Sea surface temperature [K]
rad = 'cam'; % Radiation scheme

% Figure parameters
ANOT = {'(a)','(b)'}; % Panel markers
cmap = [117 112 179; 217 150 0; 27 158 119; 231 41 138; 0 0 0]/255; % Color scale
DX = 0.185; DY = 0.75; % Subplot text annotations
f_fig = {'lw','sw','sef','adv','dmsedt'}; Nf_fig = numel(f_fig); % Fields plotted in figure
fz = 12; % Fontsize for figures
LEGLAB(2,:) = {'$\mathrm{Longwave\ }\left\langle \dot{\varphi}_{\mathrm{lw}}\right\rangle /\left\langle \varphi\right\rangle$',...
    '$\mathrm{Shortwave\ }\left\langle \dot{\varphi}_{\mathrm{sw}}\right\rangle /\left\langle \varphi\right\rangle$',...
    '$\mathrm{Surface\ flux\ }\left\langle \dot{\varphi}_{\mathrm{sf}}\right\rangle /\left\langle \varphi\right\rangle $',...
    '$\mathrm{Advection\ }\left\langle \dot{\varphi}_{\mathrm{adv}}\right\rangle /\left\langle \varphi\right\rangle $',...
    '$\mathrm{Tendency\ }\partial_{t}\left\langle \varphi\right\rangle /\left\langle \varphi\right\rangle $'};
LEGLAB(1,:) = {'$\mathrm{Longwave}\ \dot{L_{\mathrm{lw}}}$',...
    '$\mathrm{Shortwave}\ \dot{L_{\mathrm{sw}}}$',...
    '$\mathrm{Surface\ flux\ }\dot{L_{\mathrm{sf}}}$',...
    '$\mathrm{Advection}\ \dot{L_{\mathrm{adv}}}$',...
    '$\mathrm{Tendency}\ \partial_{t}L$'}; % Length scale legend
lw = 1; % Thin linewidth to better represent time variability
XLIM = [0 75]; YLIM(2,:) = [-0.49 0.49]; % x bounds and y aggregation bounds
YLAB = {'$\dot{L_{i}}\ \left[\mathrm{m/s}\right] $',['$\left\langle \',...
    'dot{\varphi}_{i}\right\rangle /\left\langle \varphi\right',...
    '\rangle \ \left[\mathrm{day^{-1}}\right] $']}; % y labels
W = 0.29; % Legend's width
if strcmp(Lmet,'Km1')==1, YLIM(1,:) = [-2.2 3]; % y length scale bounds
elseif strcmp(Lmet,'LAM')==1, YLIM(1,:) = [-7.5 7.5]; % y length scale bounds
end

load(['MAT_DATA',filesep,rad,num2str(SST),'.mat']); % Load data

%% 2. Figure
figure; set(gcf,'Position',[50 50 800 600]); % Figure's position

for ifig = 1:2, S(ifig) = subplot(2,1,ifig); % Subplot
    
    % Main plot
    line(XLIM,[0 0],'Color',[0.6 0.6 0.6],'Linewidth',lw); hold on; % Zero line for reference
    for i = 1:Nf_fig
        if ifig==1, TOPLOT = dat.(Lmet).(f_fig{i})(i2:end); % L(t)
        elseif ifig==2, TOPLOT = spd*dat.AGG.(f_fig{i})(i2:end); % AGG(t)
        end; P(i) = plot(dat.t(i2:end),movmean(squeeze(TOPLOT),WEEK),...
            'color',cmap(i,:),'Linewidth',lw); hold on; % Line plot for each field
    end; grid on;
    
    % Axes
    if ifig==2, xlabel('$\mathrm{Time\left[day\right]}$',...
            'fontsize',fz+1,'Interpreter','Latex');
    elseif ifig==1, set(gca,'XTickLabel','');
    end
    ylabel(YLAB{ifig},'Fontsize',fz+1,'Interpreter','Latex');
    set(gca,'XminorTick','on','YminorTick','on','TickLabelInterpreter','Latex');
    set(gca,'Fontsize',fz,'TitleFontWeight','normal');
    box on; drawnow; xlim(XLIM); ylim(YLIM(ifig,:));
    
    % Change subplot's and legend's positions
    if ifig==1, LP = S(ifig).Position; % Save initial S(1) position
    elseif ifig==2, S(2).Position(2) = S(1).Position(2)-S(2).Position(4);
    end
    LEG = legend(P,LEGLAB(ifig,:),'Location','eastoutside',...
        'Fontsize',fz,'Interpreter','Latex'); % Legend
    LEG.Position(2)=S(ifig).Position(2);LEG.Position(4)=S(ifig).Position(4);
    LEG.Position(1) = LP(1)+LP(3)-W; LEG.Position(3) = W;
    S(ifig).Position(3) = LEG.Position(1)-S(ifig).Position(1);
    
end

% Panel markers (a) & (b)
for i = 1:2, T(i) = annotation('textbox','String',ANOT{i},...
        'Interpreter','Latex','Color','k','Fontsize',fz+2,'Edgecolor','none');
    T(i).Position(1) = S(i).Position(1)-DX*S(i).Position(3);
    T(i).Position(2) = S(i).Position(2)+DY*S(i).Position(4);
end

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'Figure05';
elseif strcmp(Lmet,'Km1')==1, name = 'FigureS5';
end; gcfsavepdf([basedir,'PDF_DATA',filesep,name,'.pdf']); % Save plot