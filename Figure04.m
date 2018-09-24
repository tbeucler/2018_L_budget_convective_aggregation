%% Figure04.m
% tbeucler - 9/22/2018
% tbeucler - 8/30/2018
% (a) Length scale [1000 km] versus time [day] for all 20 simulations
% (b) Spatial standard deviation of MSE [MJ/m2] versus time [day]

close all; fclose('all'); clearvars;

%% 1. Parameters
% User's choice
Lmet = 'LAM'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Physical parameters
DAY = 24; % Number of timesteps in 1 day
WEEK = 7*DAY; % Number of timesteps in 1 week

% Data parameters
i31 = 24*30+1; % Beginning of thirty-first day of simulation
i75 = 24*75; % End of the seventy-fifth day of simulation
dom = 1:i75; % Total simulation time
mat = i31:i75; % Last months of simulation data
rad_array = {'rrtm','cam',''}; Nrad = numel(rad_array); % Radiation schemes
SST_array = 280:5:310; NSST = numel(SST_array); % Sea surface temperatures [K]

% Figure parameters
ANOX = 0.11; ANOY = 0.06;
cmap = parula(7); % Parula color scale (parallel with Fig6 of Wing & Cronin)
cmap(7,:) = [139 69 19]/255; % Make the yellow saddlebrown
cmapMD = [[153 51 255]/255;[153 0 153]/255;[1 0 0];[0.5 0.5 0.5];...
    [0.5 0.5 0.5];[0 0 0];[0 0 0]]; % Color scale for mechanism den. exp.
CMAP = cat(3,cmap,cmap,cmapMD); % Full color scale
FOR = {'%02.1f','%02.0f'}; % Format of each legend's line
fz = 12; % Fontsize
lw = 1; % Linewidth
STR_array = {'$\mathrm{MD_{SFC}^{RAD}}','$\mathrm{MD_{RAD}}',...
    '$\mathrm{MD_{SFC}}','$\mathrm{SQ}','','$\mathrm{BSQ}',''}; % Legends for MD lines
Xpos = [0.405 0.805]; % x-position of legends
YLAB = {'L [1000km]',['$\mathrm{Standard\ deviation\ FMSE\ \left',...
    '[MJ.m^{-2}\right]}$']}; % y-axis labels
YLOC = {'left','right'}; % Location of y-axis

% Legend for each figure line
a = zeros(2,1); b = a; YLIM = zeros(2,2); adj = zeros(2,Nrad,NSST); % Position of legends in Figure
a(2) = 1/2.5; b(2) = -1.115+1.2; % Linear position of agg. legends (a*AGG+b)
adj(2,:,:) = [[0 0 0 0 0 0 0];[-2 -3 -1 -3 -1 -2 -1];[1 -3 -1 1 0 2 0]]/1e2; % Adjustment for y-position of AGG leg.
XLIM = [0 97.5]; YLIM(2,:) = [0 1.95]; % x and y bounds of aggregation plot
Xmargin = 0.1; Ymargin = 0.1; % x and y margins of subplots
if strcmp(Lmet,'Km1')==1, a(1) = 0.18075; b(1) = 0.0682; % Leg pos (a*L+b)
    adj(1,:,:) = [[0 0 0 0 0 0 0];[-2 6 -4 -5 0 -2 0];[2 -2 0 -1 0 2 0]]/1e2; % Adjustment for y-position of L legends
    YLIM(1,:) = [0 4.53]; % y bounds for L subplot
elseif strcmp(Lmet,'LAM')==1, a(1) = 0.1; b(1) = 0.0682; % Leg pos (a*L+b)
    adj(1,:,:) = [[0 0 0 0 0 0 0];[-3 8 -2 -1 1 1 -0.5];[5 -8 0 -2 0 -2 0]]/1e2; % Adjustment for y-position of L legends
    YLIM(1,:) = [0 8.13]; % y bounds for L subplot
end

%% 2. Figure
F = figure('position',[100 100 1000 750]); % Figure
for irad = 1:Nrad, rad = rad_array{irad}; % Radiation scheme
    for iSST = 1:NSST, SST = SST_array(iSST); % Sea surface temperature
        
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
        
        % Convert aggregation metric from spectral to physical space
        if irad~=3||iSST<4, CFAC = mean(squeeze(dat.VAR.mse(1,1,DAY+1:end)./...
                dat.AGG.mse(1,1,DAY+1:end)));
        elseif iSST<7, CFAC = mean(squeeze(dat.VAR.mse(DAY+1:i75))'./...
                squeeze(dat.AGG.mse(1,1,DAY+1:i75)));
        else, CFAC = 1; % Non-existent case
        end; dat.AGG.mse = CFAC*dat.AGG.mse; % Convert MSE variance
        
        for ifig = 1:2, S(ifig) = subplot(1,2,ifig); % Subplots 1(L) 2(AGG)
            
            % Plot line
            if iSST~=7||irad~=3
                if ifig==1, TOPLOT = squeeze(dat.(Lmet).mse(:,:,dom)/1e6);
                    if irad==3&&iSST>3, TOPLOT = sqrt(2)*TOPLOT; end % Scaling factor for square domain
                elseif ifig==2, TOPLOT = log10(squeeze(sqrt(dat.AGG.mse(:,:,dom))/1e6));
                end; P = plot(dat.t(dom),movmean(TOPLOT,WEEK),'color',...
                    CMAP(iSST,:,irad),'Linewidth',lw); hold on;
            end
            
            % Axis and labels
            if strcmp(rad,'rrtm')==1, set(P,'Linestyle',':'); end
            xlim(XLIM); ylim(YLIM(ifig,:)); grid on;
            xlabel('Time [day]','Interpreter','Latex');
            ylabel(YLAB{ifig},'Interpreter','Latex');
            set(gca,'TickLabelInterpreter','Latex','Box','on',...
                'TickDir','out','Fontsize',fz,'YAxisLocation',YLOC{ifig});
            
            % Legend
            if ifig==1, MAT = TOPLOT(mat); elseif ifig==2, MAT = 10.^TOPLOT(mat); end
            if irad==2||irad==3
                if irad==2, STR = ['LC',num2str(SST),'(',num2str(nanmean(...
                        MAT),FOR{ifig}),')'];
                else, STR = [STR_array{iSST},'\left(',...
                        num2str(nanmean(MAT),FOR{ifig}),'\right)$'];
                    if iSST==5||iSST==7, STR = ''; end
                end
                Ypos = double(a(ifig)*nanmean(TOPLOT(i75-WEEK:i75))+...
                    b(ifig))+adj(ifig,irad,iSST);
                T = annotation('textbox',[Xpos(ifig) Ypos 0.11 0.06],'String',STR,...
                    'LineStyle','none','Interpreter','latex','Color',CMAP(iSST,:,irad),...
                    'FitBoxToText','off','BackgroundColor','none','Fontsize',fz);
            end
        end
    end
end

% Take advantage of the entire figure space
set(S(1),'Position',[Xmargin Ymargin 0.5-Xmargin (1-2*Ymargin)]);
set(S(2),'Position',[0.5 Ymargin 0.5-Xmargin (1-2*Ymargin)]);

% Change aggregation y-scale to logarithmic
YTIK = S(2).YTickLabel;
for iy = 1:numel(YTIK), YTIK{iy}=num2str(10^str2double(YTIK{iy}),FOR{1}); end
set(S(2),'YTickLabel',YTIK);

% Add subplot's ID letters
ANOT = {'a)','b)'};
for i = 1:2, T = annotation('textbox',[S(i).Position(1) S(i).Position(2)+...
        S(i).Position(4)-ANOY ANOX ANOY],'String',ANOT{i},'LineStyle',...
        'none','Interpreter','latex','FitBoxToText','off',...
        'BackgroundColor','none','Fontsize',1.5*fz);
end

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'Figure04';
elseif strcmp(Lmet,'Km1')==1, name = 'FigureS4';
end; gcfsavepdf([basedir,'PDF_DATA',filesep,name,'.pdf']); % Save plot