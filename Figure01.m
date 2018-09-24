%% Figure01.m
% tbeucler - 9/22/2018
% tbeucler - 8/29/2018
% tbeucler - 2/26/2018
% Hovmuller diagram of MSE anomaly field of 
% (a) MDRAD (b) MDSFC (c) MDRAD+SFC (d) LC300CAM

close all; fclose('all'); clearvars;

% User's choice
Lmet = 'LAM'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Data parameters
sim_array = {'mdrad','mdsfc','mdradsfc','cam300'}; % Simulations to plot

% Figure parameters
fz = 12; lw = 2; % Fontsize and Linewidth
Cmin = -30; Cmax = 30; % MSE spatial anomaly bounds [MJ/m2]
CLAB = ['$\mathrm{FMSE\ spatial\ anomaly}\ \left',...
    '[\mathrm{MJ\ m^{-2}}\right]$']; % Colorbar label
title_list = {'$\mathrm{\left(a\right)\ MD\ RAD} $';...
    '$\mathrm{\left(b\right)\ MD\ SFC} $';...
    '$\mathrm{\left(c\right)\ MD\ RAD+SFC} $';...
    '$\mathrm{\left(d\right)\ LC300\ CAM} $'}; % Title for each sub-plot
Xfig = 750; Yfig = 4/3*650; % Standard figure dimensions
ColorbarH = 0.025; % Colorbar vertical extent

figure('position',[50 50 Xfig Yfig]) % Figure
for ifig = 1:4, S(ifig) = subplot(4,3,3*(ifig-1)+(1:3)); % Subplot
    
    load(['MAT_DATA',filesep,sim_array{ifig},'.mat']); % Load sim. data
    
    MSEanomaly = (dat.MSE-repmat(mean(dat.MSE(:)),numel(dat.x),1))/1e6; % MJ/m2
    
    % Hovmuller plot of MSE spatial anomaly
    P = pcolor(dat.x/1e6,dat.t,MSEanomaly'); P.LineStyle = 'none';
    colormap(redblue); caxis([Cmin Cmax]); % Colorbar
    
    % Axes and labels
    ylabel('Time [day]','Fontsize',fz,'Interpreter','Latex');
    set(gca,'Fontsize',fz,'TickLabelInterpreter','Latex','Layer','top');
    if ifig==4, xlabel('x [1000km]','Fontsize',fz,'Interpreter','Latex');
    else, set(gca,'XTickLabel',''); end
    title(title_list{ifig},'Fontsize',fz,'Interpreter','Latex');
    
    % Colorbar label
    if ifig==1, c = colorbar('Northoutside');
        set(c.Label,'String',CLAB,'FontSize',fz,'Interpreter','Latex',...
            'Position',[0 2.25 0]);
    end; hold on;
    
    % Plot L(t)
    plot(squeeze(dat.(Lmet).mse)/1e6,dat.t,'color',[0 0 0],'Linewidth',lw);
    
end; hold on; 

% Subplot position's adjustment
S(4).Position(4) = S(3).Position(4);
for ifig = 1:4, S(ifig).Position(2) = 0.9-0.2*ifig; end
c.Position = [S(1).Position(1) 0.9 S(1).Position(3) ColorbarH];
c.TickLabelInterpreter = 'Latex';
S(1).Position(4) = S(2).Position(4);

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'Figure01';
elseif strcmp(Lmet,'Km1')==1, name = 'FigureS1';
end; gcfsavepdf([basedir,'PDF_DATA',filesep,name,'.pdf']); % Save plot
