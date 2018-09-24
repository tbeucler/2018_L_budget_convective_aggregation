%% FigA1.m
% tbeucler - 9/23/2018
% tbeucler - 9/18/2017
% Compares different definitions of length scale for the LC CAM 300 sim.

close all; fclose('all'); clearvars;

% Data parameters
i75 = 24*75; % End of the seventy-fifth day of simulation
rad = 'cam'; % Radiation scheme
SST = 300; % Sea surface temperature [K]

% Physical parameters
DAY = 24; % Number of timesteps in 1 day
WEEK = 7*DAY; % Number of timesteps in 1 week

% Figure parameters
a = 0.25/1e6; b = 0.0682; % Linear relation: y-pos. of legends = a*L+b
adjL = [6 -3 0 0]/1e2; % Adjustment for Y-position of L legends
COL = [[0 0 0];[125 125 255]/255;[0 0 1];[1 0 0]]; % Color of each line
F = {'LAM','Km1','GEOM','AUT'}; NF = 4; % Different definitions of L
fz = 12; % Fontsize
LAB = {['$\frac{\left\langle \lambda\varphi\',...
    'right\rangle }{\left\langle \varphi\right\rangle }$'],...
    ['$\frac{\left\langle \varphi\right\',...
    'rangle }{\left\langle k\varphi\right\rangle }$'],...
    ['$\sqrt{\frac{\left\langle \lambda\',...
    'varphi\right\rangle }{\left\langle k\varphi\right\rangle }}$'],...
    '$\ L_{\varphi}$'}; % Label for each line
lw = 2.5; % Linewidth
XposL = 0.7; % Position of legend for each line
XLIM = [0 100]; YLIM = [0 3.5]; % X and Y bounds for plotting

load(['MAT_DATA',filesep,rad,num2str(SST)]); % Load data
figure; set(gcf,'Position',[50 50 500 500]); % Figure

for iF = 1:NF, plot(dat.t(1:i75),movmean(squeeze(dat.(F{iF}).mse(1:i75)/1e6),WEEK),...
        'color',COL(iF,:),'Linewidth',lw); hold on; % Plot one line per L(t)
    % Legend for each line
    Ypos = double(a*nanmean(dat.(F{iF}).mse(i75-WEEK:i75))+b)+adjL(iF);
    T = annotation('textbox',[XposL Ypos 0.11 0.06],'String',LAB{iF},...
        'LineStyle','none','Interpreter','latex','Color',COL(iF,:),...
        'FitBoxToText','off','BackgroundColor','none','Fontsize',2*fz);
end; grid on; xlim(XLIM); ylim(YLIM); % x and y bounds for plotting
xlabel('Time [day]','Fontsize',fz,'Interpreter','Latex'); % x label
ylabel('Length scale [1000 km]','Fontsize',fz,'Interpreter','Latex');
set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex'); % Axes appearance

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir,'PDF_DATA',filesep,'FigureA1.pdf']); % Save plot