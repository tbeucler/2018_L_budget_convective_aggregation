%% Figure06.m
% tbeucler - 9/23/2018
% tbeucler - 8/17/2018
% Summary of all 20 experiments in (L rate,AGG rate) space
% Creates a summary table in parallel in the local folder TEX_DATA

close all; clearvars; fclose('all');

%% 1. Parameters and initialization of figure and table
% User's choice
Lmet = 'LAM'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Physical parameters
DAY = 24; % Number of timesteps in 1 day
spd = 24*3600; % Number of seconds per day

% Model parameters
i2 = 25; % Beginning of second day of simulation
i31 = 24*30+1; % Beginning of thirty-first day of simulation
i75 = 24*75; % End of the seventy-fifth day of simulation
gro = i2:(i31-1); mat = i31:i75; % First month and rest of the simulation
SST_array = 280:5:310; NSST = numel(SST_array); % Sea surface temperatures [K]
rad_array = {'cam','rrtm',''}; Nrad = numel(rad_array); % Radiation scheme
if strcmp(Lmet,'LAM')==1, sq_fac = 5; % Mult. factor for MD SFC/SQ/BSQ sim.
elseif strcmp(Lmet,'Km1')==1, sq_fac = 1; % Mult. factor for MD SFC/SQ/BSQ sim.
end
Xdom = 1.2288e7/2; % Half of domain size (max. resolvable scale) [m]

% Figure parameters
ARAD = {'CAM','RRTM','MD'}; % Array legend prefixes
AMD = {'MD RAD SFC','MD RAD',[num2str(sq_fac),' MD SFC'],...
    [num2str(sq_fac),' SQ CAM'],[num2str(sq_fac),' SQ RRTM'],...
    [num2str(sq_fac),'BSQ RRTM'],''}; % Array legend prefixes (continued)
colors_rgb = [117 112 179; 217 150 0; 27 158 119; 231 41 138]; % Longwave/Shortwave/SEF/Advection
cmap = colors_rgb./255; % Adapted colorbar from Wing & Cronin (2016)
delgraph = [0 3.5e-3]; % Space between figures
Dellef = 0.07; % Size for axis on the left of the graph
Delbot = 0.1; % Size for axis on the bottom of the graph
DelSST = (1-Dellef)/NSST; Delrad = (1-Delbot)/Nrad; % Size of each subplot
DelX = DelSST/2; DelY = Delrad/7.5; % Annotation's box size
dlineSST = Dellef; dlinerad = Delbot; % Line extending from box plots
F = {'Xgro','Xmat','Ygro','Ymat'}; NF = numel(F); % Plotted fields
f_fig = {'lw','sw','sef','adv'}; Nf_fig = numel(f_fig); % Fields plotted in figure
fm = '%01.1f'; % Table format for aggregation and length scale rates
fz = 15; % Fontsize
LEG_TEXT = {'Longwave - Growth ','Shortwave - Growth ',...
    'S.E.F. - Growth ','Advection - Growth ',...
    'Longwave - Mature ','Shortwave - Mature ',...
    'S.E.F. - Mature ','Advection - Mature '}; % Legend text
lw = 2.5; % Linewidth
sz = 33; % Marker size
if strcmp(Lmet,'Km1')==1, YLIM = [1.15 1.15 1.15]; % Y range of figure
    XLIM = [0.79 0.79 0.79 0.29 0.29 0.29 0.29]; % X range of figure (function of SST)
elseif strcmp(Lmet,'LAM')==1, YLIM = [1.15 1.15 1.15]; % Y range of figure
    XLIM = [0.275 0.275 0.275 0.275 0.275 0.275 0.275]; % X range of figure (function of SST)
end
XLAB = '$\dot{L}_{i}/L\ \left[\mathrm{day^{-1}}\right]$'; % X label
YLAB = ['$\left\langle \dot{\varphi}_{i}\right\rangle /',...
                '\left\langle \varphi\right\rangle \ \left[',...
                '\mathrm{day^{-1}}\right]$']; % Y label

% Create .tex file for Table 2
if ~exist('TEX_DATA','dir'), mkdir('TEX_DATA'); end % Create TEX_DATA folder if doesn't exist
filename = ['TEX_DATA',filesep,'Table_',Lmet,'.tex']; fileid = fopen(filename,'w');
% Define colors for table and set up horizontal margins to 1.5pt
for icol = 1:size(cmap,1), CSTR = num2str(cmap(icol,1:3),'%f,');
    fprintf(fileid,['\\definecolor{C',num2str(icol),'}{rgb}{',...
        CSTR(1:end-1),'}\n']);
end; fprintf(fileid,'\\setlength\\tabcolsep{1.5pt}\n\n');
% Begin table
fprintf(fileid,'\\begin{table}[H]\n'); fprintf(fileid,'{\\footnotesize\n');
fprintf(fileid,'\\begin{centering}\n');

figure
set(gcf,'Position',[10 10 1200 600]);

for irad = 1:Nrad, rad = rad_array{irad};
    for iSST = 1:NSST, SST = SST_array(iSST);
        
        % Name and load each experiment
        if irad==3
            if iSST==1, name = 'mdradsfc'; % MD RAD+SFC experiment
            elseif iSST==2, name = 'mdrad'; % MD RAD experiment
            elseif iSST==3, name = 'mdsfc'; % MD SFC experiment
            elseif iSST==4, name = 'sqcam'; % CAM square experiment
            elseif iSST==5, name = 'sqrrtm'; % RRTM square experiment
            elseif iSST==6, name = 'bsqrrtm'; % Big square experiment
            end
        else, name = [rad,num2str(SST)];
        end; load(['MAT_DATA',filesep,name,'.mat']); % Load data
        
        %% 2. Prepare reduced data for figure
        
        toplot = struct; % Structure containing time-averaged fields to plot
        for i = 1:Nf_fig
            % Aggregation tendencies averaged over growth and mature phases
            toplot.Ygro.(f_fig{i}) = spd*squeeze(mean(dat.AGG.(f_fig{i})(:,:,gro),3));
            toplot.Ymat.(f_fig{i}) = spd*squeeze(mean(dat.AGG.(f_fig{i})(:,:,mat),3));
            % Length scale tendencies averaged over growth and mature phases
            toplot.Xgro.(f_fig{i}) = spd*squeeze(mean(dat.(Lmet).(f_fig{i})...
                (:,:,gro)./dat.(Lmet).mse(:,:,gro),3));
            toplot.Xmat.(f_fig{i}) = spd*squeeze(mean(dat.(Lmet).(f_fig{i})...
                (:,:,mat)./dat.(Lmet).mse(:,:,mat),3));
            % Diabatic length scale averaged over growth and mature phases
            for iPHA = 1:2
                if iPHA==1, dom = gro; PHA='LDIABgrowth';
                elseif iPHA==2, dom = mat; PHA='LDIABmature';
                end
                if strcmp(Lmet,'Km1')==1
                    L_DIAB = nanmean(dat.(Lmet).mse(:,:,dom),3)./(1-...
                        nanmean(dat.(Lmet).(f_fig{i})(:,:,dom)./...
                        (dat.(Lmet).mse(:,:,dom)),3)./...
                        nanmean(dat.AGG.(f_fig{i})(:,:,dom),3));
                    L_DIAB(L_DIAB>Xdom) = Xdom; % Diabatic length scale cannot be larger than Xdom
                    L_DIAB(L_DIAB<0) = Xdom; % Set super L rate to max value
                elseif strcmp(Lmet,'LAM')==1
                    L_DIAB = nanmean(dat.(Lmet).mse(:,:,dom),3)+...
                        nanmean(dat.(Lmet).(f_fig{i})(:,:,dom)./...
                        nanmean(dat.AGG.(f_fig{i})(:,:,dom),3));
                end
                toplot.(PHA).(f_fig{i}) = L_DIAB/1e6; % Convert to Mm
            end
        end
        % Length scale averaged over growth and mature phases
        toplot.Lgrowth = squeeze(nanmean(dat.(Lmet).mse(:,:,gro),3))/1e6; % In Mm
        toplot.Lmature = squeeze(nanmean(dat.(Lmet).mse(:,:,mat),3))/1e6; % In Mm
        
        %% 3. Figure
        
        S(irad,iSST)=subplot(Nrad,NSST,NSST*(irad-1)+iSST);
        
        % Background quadrant
        line([-XLIM(iSST) XLIM(iSST)],[0 0],'Linewidth',lw,...
            'color',[0.6 0.6 0.6]); hold on;
        line([0 0],[-YLIM(irad) YLIM(irad)],'Linewidth',lw,...
            'color',[0.6 0.6 0.6]); hold on;
        % Plot lines from origin
        for i = 1:Nf_fig
            XG = toplot.Xgro.(f_fig{i}); YG = toplot.Ygro.(f_fig{i});
            XM = toplot.Xmat.(f_fig{i}); YM = toplot.Ymat.(f_fig{i});
            if irad==3&&iSST>2 % Multiply SQ & MD SFC feedbacks by sq_fac
                XG=sq_fac*XG; XM=sq_fac*XM; YG=sq_fac*YG;YM=sq_fac*YM;
            end
            line([0 XG],[0 YG],'Linewidth',lw,'color',cmap(i,:),...
                'Linestyle',':'); hold on;
            line([0 XM],[0 YM],'Linewidth',lw,'color',cmap(i,:)); hold on;
        end
        % Scatter forcings in (L,AGG) space
        for i = 1:Nf_fig
            XG = toplot.Xgro.(f_fig{i}); YG = toplot.Ygro.(f_fig{i});
            XM = toplot.Xmat.(f_fig{i}); YM = toplot.Ymat.(f_fig{i});
            if irad==3&&iSST>2==1 % Multiply SQ & MD SFC feedbacks by sq_fac
                XG=sq_fac*XG; XM=sq_fac*XM; YG=sq_fac*YG;YM=sq_fac*YM;
            end
            G(i)=scatter(XG,YG,'MarkerFaceColor','w','MarkerEdgeColor',...
                cmap(i,:),'SizeData',sz,'Linewidth',lw); hold on;
            M(i)=scatter(XM,YM,'MarkerFaceColor',cmap(i,:),...
                'MarkerEdgeColor',cmap(i,:),'SizeData',sz); hold on;
        end
        % Figure properties
        grid on; xlim([-XLIM(iSST) XLIM(iSST)]); ylim([-YLIM(irad) YLIM(irad)]);
        set(findobj(gcf,'type','axes'),'TickLabelInterpreter','Latex',...
            'Fontsize',fz,'Linewidth',lw);
        if iSST>1, set(gca,'YTickLabel',''); end
        if irad<3, set(gca,'XTickLabel',''); end
        if iSST==1, ylabel(YLAB,'Fontsize',fz,'Interpreter','Latex'); end
        if irad==3, xlabel(XLAB,'Fontsize',fz,'Interpreter','Latex'); end
        TAB.(['SST',num2str(iSST)]) = toplot; % Save structure (SST) for Table
        
    end
    
    %% 4. Table
    if irad==3, NSST = NSST-1; end % One less MD experiments than LC
    % 4.1 Print top line of table
    fprintf(fileid,'\\begin{tabular}{|c|c|c|c|c|c|c|c|}\n');
    fprintf(fileid,'\\hline\n');
    fprintf(fileid,'Var$\\downarrow\\ $Exp$\\rightarrow$ ');
    for iSST = 1:NSST, SST = SST_array(iSST);
        if irad~=3, fprintf(fileid,['& ',ARAD{irad},num2str(SST),' ']);
        else, fprintf(fileid,['& ',AMD{iSST},' ']);
        end
    end; fprintf(fileid,'\\tabularnewline\n');
    fprintf(fileid,'\\hline\n'); fprintf(fileid,'\\hline\n');
    
    % 4.2 Growth and mature phase's length scale
    fprintf(fileid,'$L_{\\mathrm{G/M}}\\ \\left[1000\\mathrm{km}\\right]$ ');
    for iSST = 1:NSST
        fprintf(fileid,['& ',num2str(TAB.(['SST',num2str(iSST)]).Lgrowth,fm),...
            '/',num2str(TAB.(['SST',num2str(iSST)]).Lmature,fm),'\\ ']);
    end; fprintf(fileid,'\\tabularnewline\n'); fprintf(fileid,'\\hline\n');
    
    % 4.3 Growth and mature phase's aggregation rate
    for iphase = 1:2
        if iphase==1, PHA='G'; else, PHA='M'; end
        fprintf(fileid,['$\\left\\langle \\dot{\\varphi}_{i}\\right\\rangle /',...
            '\\left\\langle \\varphi\\right\\rangle _{\\mathrm{',PHA,...
            '}}\\ \\left[\\mathrm{d^{-1}}\\right]$ ']);
        for iSST = 1:NSST, fprintf(fileid,'& ');
            for i = 1:Nf_fig
                fprintf(fileid,['\\textcolor{C',num2str(i),'}{',num2str(...
                    TAB.(['SST',num2str(iSST)]).(F{iphase+2}).(f_fig{i}),fm),'}\\ ']);
            end; fprintf(fileid,' ');
        end; fprintf(fileid,'\\tabularnewline\n'); fprintf(fileid,'\\hline\n');
    end
    
    % 4.4 Growth and mature phase's length scale rate
    for iphase = 1:2
        if iphase==1, PHA='G'; else, PHA='M'; end
        fprintf(fileid,['$10\\ \\dot{L}_{i}/L_{\\mathrm{',PHA,'}}\\ \\left[',...
            '\\mathrm{d^{-1}}\\right]$ ']);
        for iSST = 1:NSST, fprintf(fileid,'& ');
            for i = 1:Nf_fig
                fprintf(fileid,['\\textcolor{C',num2str(i),'}{',num2str(...
                    sq_fac*TAB.(['SST',num2str(iSST)]).(F{iphase}).(f_fig{i}),fm),'}\\ ']);
            end; fprintf(fileid,' ');
        end; fprintf(fileid,'\\tabularnewline\n'); fprintf(fileid,'\\hline\n');
    end
    
    % 4.5 Growth and mature phase's diabatic length scale
    for iphase = 1:2
        if iphase==1, PHA='G'; STR = 'LDIABgrowth';
        else, PHA='M'; STR = 'LDIABmature';
        end
        fprintf(fileid,['${\\cal L}_{i,\\mathrm{',PHA,'}}\\ \\',...
            'left[1000\\mathrm{km}\\right]$ ']);
        for iSST = 1:NSST, fprintf(fileid,'& ');
            for i = 1:Nf_fig
                fprintf(fileid,['\\textcolor{C',num2str(i),'}{',num2str(...
                    TAB.(['SST',num2str(iSST)]).(STR).(f_fig{i}),fm),'}\\ ']);
            end; fprintf(fileid,' ');
        end; fprintf(fileid,'\\tabularnewline\n'); fprintf(fileid,'\\hline\n');
    end
    
    % 4.6 Finish tabular
    if irad~=3, fprintf(fileid,' &  &  &  &  &  &  & \\tabularnewline\n');
    else,  fprintf(fileid,' &  &  &  &  &  & \\tabularnewline\n');
    end
    fprintf(fileid,'\\hline\n');
    fprintf(fileid,'\\end{tabular}\n');
    
end; LEG = legend([G M],LEG_TEXT); NSST = NSST+1; % Replace NSST w/ or. value

% Finish table
fprintf(fileid,'\\par\\end{centering} }\n\\end{table}\n');

% Adjust figure's subplots

for irad = 1:Nrad
    for iSST = 1:NSST
        S(irad,iSST).Position(1) = Dellef+(iSST-1)*DelSST+delgraph(1);
        S(irad,iSST).Position(2) = Delbot+(Nrad-irad)*Delrad+delgraph(2);
        S(irad,iSST).Position(3) = DelSST-delgraph(1);
        S(irad,iSST).Position(4) = Delrad-delgraph(2);
        if irad<=2, ANOT = [ARAD{irad},' ',num2str(SST_array(iSST))];
        else, ANOT = AMD{iSST}; end
        if irad==3&&iSST==1,DY=3*DelY; else, DY=2*DelY; end
        if irad==3&&iSST>2, AX=S(irad,iSST).Position(1);
        else AX=S(irad,iSST).Position(1)+S(irad,iSST).Position(3)-DelX;
        end
        if irad~=Nrad || iSST~=NSST, A(irad,iSST) = annotation('textbox',...
                [AX S(irad,iSST).Position(2) DelX DY],'Fontsize',fz,'String',ANOT,...
                'BackgroundColor',[1 1 1],'Interpreter','Latex','Linewidth',lw,...
                'HorizontalAlignment','center'); hold on;
        end
    end
end
set(LEG,'Position',S(Nrad,NSST).Position,'Fontsize',fz/1.6,'Interpreter','Latex');

for irad = 1:Nrad
    L(irad) = annotation('line',...
        'X',[S(irad,1).Position(1)-dlineSST S(irad,1).Position(1)],...
        'Y',[S(irad,1).Position(2) S(irad,1).Position(2)],'Linewidth',lw,...
        'color','k');
end
for iSST = 1:NSST
    L(iSST) = annotation('line',...
        'X',[S(3,iSST).Position(1) S(3,iSST).Position(1)],...
        'Y',[S(3,iSST).Position(2)-dlinerad S(3,iSST).Position(2)],...
        'Linewidth',lw,'color','k');
end

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'Figure06';
elseif strcmp(Lmet,'Km1')==1, name = 'FigureS6';
end; gcfsavepdf([basedir,'PDF_DATA',filesep,name,'.pdf']); % Save plot
