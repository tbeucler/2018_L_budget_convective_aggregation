%% Figure03.m
% tbeucler - 9/22/2018
% Plots the octave-averaged MSE power spectrum for the reference
% Normalized by the wavenumber-averaged MSE power spectrum
% LC300CAM simulation, plotted versus wavelength
% Time-averaged over the first 30 days

fclose('all'); close all; clearvars;

%% 1. Parameters
% User's choice
Lmet = 'Km1'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Physical parameters
DAY = 24; % Number of timesteps in 1 day
WEEK = 7*DAY; % Number of timesteps in 1 week
spd = 24*3600; % Number of seconds per day

% Model parameters
rad = 'cam'; % Radiation scheme
SST = 300; % Sea surface temperature

% Data parameters
i2 = 25; % Beginning of second day of simulation
i31 = 24*30+1; % Beginning of thirty-first day of simulation 
gro = i2:i31; % First month of simulation data

% Figure parameters
colLW = [117 112 179]/255; % Color used for longwave flux in paper
fz = 12; % Fontsize for figures (default = 5)
hl = 5; % Headlength of arrows
LEGX = 0.275; LEGDX = 0.2; LEGDY = 0.075; % Legend width & height
LEGY = [0.85 0.8 0.65 0.6 0.45 0.4]; % Legend's y positions
lw = 2.5; % Linewidth
sz = 50; % Marker size
Xdom = 1.2288e7/2; % Half of domain size (max. resolvable scale) [m]
XLIM = [3.5 7]; YLIM = [-12 50]; % x and y limits 
Y_A = YLIM(2)-YLIM(1); DIVY = 7; % Vertical extent of L in figure
Yarrow = [0.25 0.5 0.75]; % Position of three arrows from L to L_LW

%% 2. Data preparation
load(['MAT_DATA',filesep,rad,num2str(SST),'.mat']); % Load data

% Decomposition of length scale tendency
AGG = nanmean(dat.AGG.lw(:,:,gro),3); % Time-averaged aggregation rate [1/s]
L = nanmean(dat.(Lmet).mse(:,:,gro),3); % Time-averaged length scale [m]
if strcmp(Lmet,'LAM')==1
    LLW = nanmean(dat.(Lmet).lw(:,:,gro),3); % L tendency [m/s]
    L_DIAB = L+LLW/AGG; % Diabatic length scale from agg. and L. rates [m]
elseif strcmp(Lmet,'Km1')==1
    LLW = nanmean(dat.(Lmet).lw(:,:,gro)./(dat.(Lmet).mse(:,:,gro)),3); % L rate
    L_DIAB = L/(1-LLW/AGG); % Diabatic length scale from logistic factor
    L_DIAB(L_DIAB>Xdom) = Xdom; % Diabatic length scale cannot be larger than Xdom
    L_DIAB(L_DIAB<0) = Xdom; % Set super shrinking/stretching rate to max value
end

% Latex legend for each term of the figure
LEGL = {'$L$','${\cal L}_{\mathrm{lw}}$',...
    ['$',num2str(L/1e6,'%03.2f'),'$'],...
    ['$',num2str(L_DIAB/1e6,'%03.2f'),'$']};
if strcmp(Lmet,'Km1')==1
    LEG = {['$\left(1\right)\ \mathrm{Logistic\ growth\ from\ }L\ \mathrm{to}',...
        '\ {\cal L}_{\mathrm{lw}}>L$'],['$\left(1-L/{\cal L}_{\mathrm{lw}}\right)=',...
        num2str(LLW/AGG,'%03.2f'),'$'],['$\left(2\right)\ \mathrm{Aggregation\ \left',...
        '(\lambda^{-1}-\mathrm{averaged}\right)\ rate}:$'],['$\left\langle \dot',...
        '{\varphi}_{\mathrm{lw}}\right\rangle /\left\langle \varphi\right\rangle =',...
        num2str(spd*AGG,'%03.2f'),'\mathrm{d^{-1}}$'],['$\left(3\right)=',...
        '\left(1\right)\times\left(2\right)\ \mathrm{Stretching\ rate}:$'],...
        ['$\ \dot{L_{\mathrm{lw}}}/L=',num2str(spd*LLW,'%03.2f'),'\mathrm{d^{-1}}$']};
elseif strcmp(Lmet,'LAM')==1
    LEG = {['$\left(1\right)\ \mathrm{Linear\ growth\ from\ }L\ \mathrm{to}',...
        '\ {\cal L}_{\mathrm{lw}}>L$'],['${\cal L}_{\mathrm{lw}}-L=',...
        num2str((L_DIAB-L)/1e3,'%03.0f'),'\mathrm{km}$'],...
        ['$\left(2\right)\ \mathrm{Aggregation\ \left',...
        '(\lambda^{-1}-\mathrm{averaged}\right)\ rate}:$'],['$\left\langle \dot',...
        '{\varphi}_{\mathrm{lw}}\right\rangle /\left\langle \varphi\right\rangle =',...
        num2str(spd*AGG,'%03.2f'),'\mathrm{d^{-1}}$'],['$\left(3\right)=',...
        '\left(1\right)\times\left(2\right)\ \mathrm{Stretching\ velocity}:$'],...
        ['$\ \dot{L_{\mathrm{lw}}}=',num2str(LLW,'%03.2f'),'\mathrm{m s^{-1}}$']};
end; NLEG = numel(LEG); % Number of textboxes

% Normalized spectral tendency: phi/<phi>
Iphi = trapz(1./dat.lambda,dat.PHI.mse,1)./trapz(1./dat.lambda,dat.lambda.^0,1);
LWsp = spd*nanmean(dat.PHI.lw(:,:,gro)./repmat(Iphi(1,1,gro),12,1,1),3);

%% 3. Figure
% Figure: Scatter spectral aggregation tendency
figure; set(gcf,'Position',[50 50 600 400]);
line(XLIM,[0 0],'Linewidth',lw/2,'Color',[0.6 0.6 0.6]); hold on;
scatter(log10(dat.lambda),LWsp,'MarkerFaceColor',colLW,'MarkerEdgeColor',...
    colLW,'SizeData',sz,'Linewidth',lw); hold on; ylim(YLIM); grid on;

% Axes
xlim(XLIM); G = gca; XTIK = G.XTickLabel;
for ix = 1:numel(XTIK), XTIK{ix}=strcat('$10^{',XTIK{ix},'}$'); end
set(gca,'XTickLabel',XTIK,'TickLabelInterpreter','Latex','Fontsize',fz);
xlabel('$\lambda\left[\mathrm{m}\right]$','Fontsize',fz,'Interpreter','Latex');
ylabel(['$\dot{\varphi}_{\mathrm{lw}}/\left\langle \varphi\',...
    'right\rangle \ \left[\mathrm{1/\mathrm{day}}\right]$'],...
    'Fontsize',fz+1,'Interpreter','Latex');

% Lines to represent L and L_LW
line([log10(L) log10(L)],[YLIM(1)+(DIVY-1)*Y_A/DIVY YLIM(2)],...
    'color','k','Linewidth',lw); hold on;
line([log10(L_DIAB) log10(L_DIAB)],[YLIM(1)+(DIVY-1)*Y_A/DIVY YLIM(2)],...
    'color',colLW,'Linewidth',lw); hold on;

% Arrows between the lines
Scaling_X = G.Position(3)/(G.XLim(2)-G.XLim(1));
Scaling_Y = G.Position(4)/(G.YLim(2)-G.YLim(1));
XA_left = (log10(L)-G.XLim(1))*Scaling_X + G.Position(1);
XA_right = (log10(L_DIAB)-G.XLim(1))*Scaling_X + G.Position(1);
for i=1:3, iar = Yarrow(i); A(i) = annotation('arrow','HeadLength',hl);
    YA(i) = (YLIM(2)-iar*Y_A/DIVY-G.YLim(1))*Scaling_Y+G.Position(2);
    set(A(i),'X',[XA_left XA_right],'Y',[YA(i) YA(i)],'Color',colLW);
end

% L and L_LW legend for each line
for i=1:4, LLEG(i) = annotation('textbox');
    if i==1||i==3, X=(log10(L)-0.45-G.XLim(1))*Scaling_X+G.Position(1);
    elseif i==2||i==4, X=(log10(L_DIAB)+0.1-G.XLim(1))*Scaling_X+G.Position(1);
    end; set(LLEG(i),'String',LEGL{i},'Interpreter','Latex',...
        'Fontsize',fz+2,'Edgecolor','none','Horizontalalignment','center');
    LLEG(i).Position(1) = X; LLEG(i).Position(2) = YA(1+2*floor((i-1)/2));
    LLEG(i).Position(4) = abs(YA(2)-YA(3));
    if i==1||i==3, LLEG(i).Color=[0 0 0]; else, LLEG(i).Color=colLW; end
end

% Textboxes
for i=1:NLEG, T(i) = annotation('textbox');
    set(T(i),'String',LEG{i},'Fontsize',fz,'Interpreter','Latex','color',colLW,...
        'Edgecolor','none','Position',[LEGX LEGY(i) LEGDX LEGDY],...
        'HorizontalAlignment','center');
end

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'Figure03';
elseif strcmp(Lmet,'Km1')==1, name = 'FigureS3';
end; gcfsavepdf([basedir,'PDF_DATA',filesep,name,'.pdf']); % Save plot
