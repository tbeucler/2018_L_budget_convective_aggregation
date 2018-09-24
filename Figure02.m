%% Figure02.m
% tbeucler - 9/22/2018
% tbeucler - 9/18/2018
% Color plots of MSE spatial anomalies and the corresponding length scale
% Uses the three simulations SQ CAM, SQ RRTM, and BSQ RRTM
% and 5 time-frames: 0-5,20-25,40-45,60-65,80-85 days

fclose('all'); close all; clearvars;

% User's choice
Lmet = 'LAM'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Figure parameters
Cmin = -50; Cmax = 50; % Bounds for the MSE spatial anomaly [MJ/m2]
CLAB = '$\mathrm{FMSE\ spatial\ anomaly}\ \left[\mathrm{MJ\ m^{-2}}\right]$'; % Colorbar label
dx = 0.04; dy = 0.025; % x & y distance between subplots
fz = 11; % Fontsize
lw = 1.5; % Linewidth
title_array = {'Day 0-5','Day 20-25','Day 40-45','Day 60-65','Day 80-85'};
ylab = {'SQ CAM','SQ RRTM','BSQ RRTM'}; % y labels

% Data parameters
sim_array = {'sqcam','sqrrtm','bsqrrtm'}; Nsim = numel(sim_array); % Sim. to plot
t_array = [1 5 9 13 17]; Nt = numel(t_array); % Time periods to plot

figure
set(gcf,'Position',[10 10 200*Nt 200*Nsim]); % Figure dimensions
for ifig = 1:Nsim, load(['MAT_DATA',filesep,sim_array{ifig},'.mat']); % Load data
    
    % Averages length scale over appropriate time periods
    tmin = zeros(Nt,1); tmax = tmin; % Initialization
    for it = 1:Nt % Loops over time periods to find L indices
        [~,imin]=min(abs(dat.t-5*(t_array(it)-1))); tmin(it)=imin;
        [~,imax]=min(abs(dat.t-5*(t_array(it)))); tmax(it)=imax;
    end
    
    for it = 1:Nt, S(ifig,it) = subplot(Nsim,Nt,Nt*(ifig-1)+it); % Subplots
        % Prepare MSE anomaly and L fields for plotting
        MSEanomaly = (dat.MSE(:,:,t_array(it))-repmat(mean(mean(dat.MSE(:,:,...
            t_array(it)),2),1),numel(dat.x),numel(dat.y),1))/1e6; % MJ/m2
        L = sqrt(2)*mean(dat.(Lmet).mse(:,:,tmin(it):tmax(it)),3); % m
        
        % Plot MSE anomaly field
        P = pcolor(dat.x/1e6,dat.y/1e6,MSEanomaly'); hold on;
        P.LineStyle = 'none'; colormap(redblue); caxis([Cmin Cmax]); % Colorbar
        
        % Plot L as a circle of radius L/2
        CIRCLE = circle(mean(dat.x(:))/1e6,mean(dat.y(:))/1e6,L/2e6);
        
        % Axes and labels
        if ifig==1, title(title_array{it},'Fontsize',fz,'Interpreter','Latex'); end
        if ifig<Nsim, set(gca,'XTickLabel',''); end
        if it==1, ylabel(ylab{ifig},'Fontsize',fz,'Interpreter','Latex');
        elseif it==Nt, ylabel('1000 km','Fontsize',fz,'Interpreter','Latex');
            set(gca,'YAxisLocation','right');
        else, set(gca,'YTickLabel','');
        end
        set(gca,'Fontsize',fz-1,'TickLabelInterpreter','Latex',...
            'Box','on','Boxstyle','full','color','k','Linewidth',lw);
    end
    
end

% Colobar
c = colorbar('Southoutside'); set(c.Label,'String',CLAB,'FontSize',fz,...
    'Interpreter','Latex','Position',[0 -2.25 0]);
c.Position = [0.1 0.075 0.8 0.02]; c.TickLabelInterpreter = 'Latex';

% Change subplot's positions
for ifig = 1:Nsim
    for it = 1:Nt
        S(ifig,it).Position = [(it-1)*(Nsim/(Nt*(Nsim+1)))+(it+0.2)*dx ...
            0.15+(3-ifig)*(((Nsim+1)^(-1))+dy) (Nsim/(Nt*(Nsim+1))) ((Nsim+1)^(-1))];
    end
end

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'Figure02';
elseif strcmp(Lmet,'Km1')==1, name = 'FigureS2';
end; gcfsavepdf([basedir,'PDF_DATA',filesep,name,'.pdf']); % Save plot

% Functions used within script
function h = circle(x,y,r) % Function to plot a circle
hold on; th = 0:pi/50:2*pi; % theta coordinate
xunit = r * cos(th) + x; yunit = r * sin(th) + y; % x and y coordinates
h = plot(xunit, yunit,'color',[1 1 1],'Linewidth',1.5); hold off;
end