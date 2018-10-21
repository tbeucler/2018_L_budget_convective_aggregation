%% Fig08.m
% tbeucler - 9/23/2018
% tbeucler - 8/30/2018

close all; fclose('all'); clearvars;

%% 1. Parameters
% User's choice
Lmet = 'LAM'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Physical parameters
DAY = 24; % Number of timesteps in 1 day
WEEK = 7*DAY; % Number of timesteps in 1 week
spd = 24*3600; % Number of seconds per day

% Model parameters
i75 = 24*75; % End of the seventy-fifth day of simulation
rad = 'cam'; % Radiation scheme
SST = 300; % Sea surface temperature [K]

% Figure parameters
ANOT = {'$\left(\mathrm{a}\right)$','$\left(\mathrm{b}\right)$',...
    '$\mathrm{\left(c\right)=\left(a\right)\times\left(b\right)}$'}; % Panel legends
cmap = [27 158 119; 51 51 255; 255 170 0; 230 0 255]; cmap = cmap./255; % Colorbar
FSUB = {'AGG','','Lplot'}; % Fields plotted in each subplot
f_fig = {'sef','sefh','sefu','sefe'}; Nf_fig = numel(f_fig);
fz = 12; % Fontsize for figures (default = 5)
LEGarray = {'$\mathrm{Total\ surface\ flux}$',...
    '$\mathrm{Enthalpy\ disequilibrium}$',...
    '$\mathrm{Surface\ wind\ speed}$','$\mathrm{Eddy\ term}$'}; % Legend for each subplot
lw = 2.5; % Linewidth
Ndiv = 5; % Time division for logistic term
sz = 100; % Marker size
Tx = [-0.05 -0.05 0.065]; Ty = 2.7; % Textbox reduction factors
XLAB = '$\mathrm{Time\left[day\right]}$'; % x-label
XLIM = [0 74]; % Bounds for the x axis
if strcmp(Lmet,'Km1')==1
Ymin = [-0.15 -1.1 -0.11]; Ymax = [1 0.5 0.6]; % X and Y bounds
YLAB = {['$\left\langle \dot{\varphi}_{i}\right\rangle /\left',...
    '\langle \varphi\right\rangle \ \left[\mathrm{1/day}\right]$'],...
    '$1-L/{\cal L}_{i}$','$\dot{L_{i}}/L\ \left[\mathrm{1/day}\right]$'};
elseif strcmp(Lmet,'LAM')==1
    Ymin = [-0.15 -1.6 -3.1]; Ymax = [1 0.5 5]; % X and Y bounds
YLAB = {['$\left\langle \dot{\varphi}_{i}\right\rangle /\left',...
    '\langle \varphi\right\rangle \ \left[\mathrm{1/day}\right]$'],...
    '${\cal L}_{\mathrm{i}}-L\ \left[\mathrm{1000\ km}\right]$',...
    '$\dot{L_{i}} \left[\mathrm{m s^{-1}}\right]$'}; % Y label
end

load(['MAT_DATA',filesep,rad,num2str(SST),'.mat']); % Load data

%% 2. Prepare reduced data for plotting 

% Decomposition of length scale tendency for each of the Ndiv time slices
TIME = zeros(Ndiv,1); % Initialize time coordinate
for i = 1:Nf_fig, LOG.(f_fig{i}) = zeros(Ndiv,1); % Logistic term
    for it = 1:Ndiv, dom = (1+(it-1)*floor(i75/Ndiv)):(it*floor(i75/Ndiv));
        AGG = nanmean(dat.AGG.(f_fig{i})(:,:,dom),3); % Aggregation rate
        if strcmp(Lmet,'Km1')==1, LLW = nanmean(dat.(Lmet).(f_fig{i})...
                (:,:,dom)./(dat.(Lmet).mse(:,:,dom)),3); % L rate
        L = nanmean(dat.(Lmet).mse(:,:,dom),3); % Average L
        LOG.(f_fig{i})(it) = LLW/AGG; % LW logistic term (dimensionless)
        dat.Lplot.(f_fig{i}) = spd*dat.(Lmet).(f_fig{i})./dat.(Lmet).mse; % Conversion for plotting
        elseif strcmp(Lmet,'LAM')==1
        LLW = nanmean(dat.(Lmet).(f_fig{i})(:,:,dom),3); % L rate
        LOG.(f_fig{i})(it) = (LLW/AGG)/1e6; % Linear growth rate [1000km]
        dat.Lplot.(f_fig{i}) = dat.(Lmet).(f_fig{i}); 
        end
        TIME(it) = mean(dat.t(dom)); % Time coordinate
    end; dat.AGG.(f_fig{i}) = spd*dat.AGG.(f_fig{i}); % Conversion for plotting
end

figure; set(gcf,'Position',[50 50 800 800]); % Figure

for ifig = 1:3, S(ifig) = subplot(3,1,ifig); % Subplot
    
    if ifig~=2 % Plot AGG(t) and L(t)
        for iF = 1:Nf_fig, plot(dat.t,movmean(squeeze(dat.(FSUB{ifig}).(...
                f_fig{iF})),WEEK),'Linewidth',lw,'color',cmap(iF,:)); hold on;
        end
    else % Plot arrows indicating stretching/shrinking length scale factor
        line([min(dat.t) max(dat.t)],[0 0],'color',[0.6 0.6 0.6],...
            'Linewidth',lw); hold on; % Zero line
        for iF = 3:-1:2 
            for it = 1:Ndiv
            AR = scatter(TIME(it),LOG.(f_fig{iF})(it),'MarkerEdgeColor',cmap(iF,:),...
                'MarkerFaceColor',cmap(iF,:),'SizeData',sz); hold on; % Arrow's head
            line([TIME(it) TIME(it)],[0 LOG.(f_fig{iF})(it)],...
                'Linewidth',2*lw,'Color',cmap(iF,:)); hold on; % Arrow's line
            if LOG.(f_fig{iF})(it)>0, set(AR,'Marker','^');  % Arrow up for stretching
            else, set(AR,'Marker','v'); end % Arrow down for shrinking
            end
        end
    end; grid on; xlim(XLIM); ylim([Ymin(ifig) Ymax(ifig)]); % x and y bounds
    set(gca,'XminorTick','on','YminorTick','on','TickLabelInterpreter',...
        'Latex','TitleFontWeight','normal','Fontsize',fz); box on; drawnow;
    if ifig==3, xlabel(XLAB,'fontsize',fz+1,'Interpreter','Latex'); % X label
    else, set(gca,'XTickLabel',''); end % No x label for subplots 1 and 2
    ylabel(YLAB{ifig},'Fontsize',fz+1,'Interpreter','Latex'); % Y label
    
end

LEG = legend(LEGarray,'Location','southoutside','orientation','horizontal',...
    'Fontsize',fz,'Interpreter','Latex'); % Legend

% Adjust subplot and legend positions
for ifig = 1:3, S(ifig).Position(2) = S(ifig).Position(2)-0.04*(4-ifig); end
LEG.Position(2) = S(1).Position(2)+1.1*S(1).Position(4);
LEG.Position(1) = S(1).Position(1);
LEG.Position(3) = S(1).Position(3);
S3.Position(3) = S(2).Position(3);

% Subplot annotations
for ifig = 1:3, T(ifig) = annotation('textbox','String',ANOT{ifig},...
        'Interpreter','Latex','Fontsize',fz+2,'Edgecolor','k',...
        'Color','k','BackgroundColor',[1 1 1]); % (a)/(b)/(c) subplot annotations
    T(ifig).Position(3) = T(ifig).Position(3)+Tx(ifig);
    T(ifig).Position(4) = T(ifig).Position(4)/Ty;
    T(ifig).Position(1) = S(ifig).Position(1)+S(ifig).Position(3)-T(ifig).Position(3);
    T(ifig).Position(2) = S(ifig).Position(2)+S(ifig).Position(4)-T(ifig).Position(4);
end

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'Figure08';
elseif strcmp(Lmet,'Km1')==1, name = 'FigureS8';
end; gcfsavepdf([basedir,'PDF_DATA',filesep,name,'.pdf']); % Save plot