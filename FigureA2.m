%% FigA2.m
% tbeucler - 9/23/2018
% tbeucler - 8/20/2018
% Octave-averaged normalized power spectral tendencies of
% (a) diabatic fluxes (b) ice,liquid,cloud water (c) surface enthalpy fluxes

close all; fclose('all'); clearvars;

%% 1. Parameters
% User's choice
Lmet = 'LAM'; % L metric (LAM=Phi-averaged wavelength, K^-1=Integral MSE scale)

% Physical parameters
Lv = 2.5e6; % Latent heat of vaporization of water vapor (J/kg)
spd = 24*3600; % Number of seconds per day

% Model parameters
i2 = 25; % Beginning of second day of simulation
i31 = 24*30+1; % Beginning of thirty-first day of simulation
i75 = 24*75; % End of the seventy-fifth day of simulation
gro = i2:(i31-1); mat = i31:i75; % First month and rest of the simulation
SST = 300; % Sea surface temperature [K]
rad = 'cam'; % Radiation scheme
Xdom = 1.2288e7/2; % Maximum length scale that can be resolved by set-up [m]

% Figure parameters
cmap = [117 112 179; 217 150 0; 27 158 119; 231 41 138;... % LW/SW/SEF/ADV
    153 153 153;217 150 0;117 112 179;... % QN/LIQ/ICE
    27 158 119; 51 51 255; 255 170 0; 230 0 255]/255; % Color scale for SEF
DelX = 0.0275; DelY = 0.05; % Annotation's box sizes
dxleg = 0.005; % Space for printing legend's right boundary
f_fig = {'lw','sw','sef','adv',... % Diabatic fluxes
    'qn','liq','ice',... % Fields plotted in cloud spectra
    'sef','sefh','sefu','sefe'}; Nf = numel(f_fig); % Surface enthalpy fluxes
diab = 1:4; cloud = 5:7; sef = 8:11; % Indices corresponding to each field group
LEG = {'Longwave - Growth phase','Shortwave - Growth phase',...
    'Surface enthalpy flux - Growth phase','Advection - Growth phase',...
    'Longwave - Mature phase','Shortwave - Mature phase',...
    'Surface enthalpy flux - Mature phase','Advection - Mature phase',...
    'Total cloud water - Growth phase','Liquid water - Growth phase',...
    'Ice water - Growth phase','Total cloud water - Mature phase',...
    'Liquid water - Mature phase','Ice water - Mature phase',...
    'Total surface flux - Growth phase','Enthalpy disequilibrium - Growth phase',...
    'Surface wind speed - Growth phase','Eddy term - Growth phase',...
    'Total surface flux - Mature phase','Enthalpy disequilibrium - Mature phase',...
    'Surface wind speed - Mature phase','Eddy term - Mature phase'};
xleft = 0.05; xfig = 0.575; xleg = 0.05; yfig = 0.25; ybot = (1-3*yfig)/3;
fz = 15; % Fontsize
lw = 2.5; % Linewidth
sz = 50; % Marker size
XLIM = [3.75 7]; % X bounds
XLAB = '$\lambda\left[\mathrm{m}\right]$'; % x label
YLIM = [[-0.075 0.075];[0 5e-3];[-0.035 0.055]]; % y bounds for subplots
YL = [0.75 1]; % y-line to plot the diabatic length scale [%]
YLAB = {'$\lambda^{-1}\dot{\varphi}_{i}/\left\langle \varphi\right\rangle \ \left[\mathrm{day^{-1}\ m^{-1}}\right]$',...
    '$\lambda^{-1}L_{v}\mathrm{Re}\left(\widehat{H}^{*}\widehat{q_{n}}\right)/\left\langle \varphi\right\rangle$',...
    '$\lambda^{-1}\dot{\varphi}_{i}/\left\langle \varphi\right\rangle \ \left[\mathrm{day^{-1}\ m^{-1}}\right]$'};

load(['MAT_DATA',filesep,rad,num2str(SST)]); Nlam = numel(dat.lambda); % Load data

%% 2. Prepare reduced data for figure

Iphi = trapz(1./dat.lambda,dat.PHI.mse,1); % Find normalization constant int{phi d(1/lambda)})

for i = 1:Nf % Loop over fields to plot
    if strcmp(f_fig{i},'qn')==1||strcmp(f_fig{i},'liq')==1||...
            strcmp(f_fig{i},'ice')==1, pre = Lv/2; % Pre-factor Lv for clouds
    else, pre = spd; % Convert from [1/s] to [1/d] for diabatic fluxes
    end
    % Spatial spectral tendencies averaged over growth and mature phases
    toplot.AggG.(f_fig{i}) = pre*nanmean(dat.PHI.(f_fig{i})(:,:,gro)./(...
        repmat(dat.lambda,1,1,numel(gro)).*repmat(Iphi(1,1,gro),Nlam,1,1).^1),3);
    toplot.AggM.(f_fig{i}) = pre*nanmean(dat.PHI.(f_fig{i})(:,:,mat)./(...
        repmat(dat.lambda,1,1,numel(mat)).*repmat(Iphi(1,1,mat),Nlam,1,1).^1),3);
    % Diabatic length scale averaged over growth and mature phases
    for iF = 1:2
        if iF==1, dom = gro; F='LG'; elseif iF==2, dom = mat; F='LM'; end
        if strcmp(Lmet,'Km1')==1
            L_DIAB = nanmean(dat.(Lmet).mse(:,:,dom),3)./(1-...
                nanmean(dat.(Lmet).(f_fig{i})(:,:,dom)./...
                (dat.(Lmet).mse(:,:,dom)),3)./...
                nanmean(dat.AGG.(f_fig{i})(:,:,dom),3)); % Diab. L from AGG and L rates
            L_DIAB(L_DIAB>Xdom) = Xdom; % Diabatic length scale cannot be larger than Xmax
            L_DIAB(L_DIAB<0) = Xdom; % Set super L rate to max value
        elseif strcmp(Lmet,'LAM')==1, L_DIAB = nanmean(dat.(Lmet).mse(:,:,dom),3)+...
                nanmean(dat.(Lmet).(f_fig{i})(:,:,dom))./...
                nanmean(dat.AGG.(f_fig{i})(:,:,dom),3); % Diab. L from AGG and L. rates
        end
        toplot.(F).(f_fig{i}) = log10(L_DIAB); % Convert to log scale for plotting
    end
end
toplot.LG.mse = log10(nanmean(dat.(Lmet).mse(:,:,gro),3)); % Convert to log sc for plot
toplot.LM.mse = log10(nanmean(dat.(Lmet).mse(:,:,mat),3)); % Convert to log sc for plot

%% 3. Figure

figure; set(gcf,'Position',[10 10 1200 700]); % Figure dimensions
for ifig = 1:3, S(ifig) = subplot(3,1,ifig); % Subplot
    
    if ifig==1, FIE=diab; elseif ifig==2, FIE=cloud; else, FIE=sef; end
    
    line([toplot.LG.mse toplot.LG.mse],YL*YLIM(ifig,2),...
        'Linewidth',1.5*lw,'Linestyle',':','color','k'); hold on; % L during first month
    line([toplot.LM.mse toplot.LM.mse],YL*YLIM(ifig,2),...
        'Linewidth',1.5*lw,'color','k'); hold on; % L during following months
    
    for i = FIE(1):FIE(end) % Loop over fields to plot
        plot(log10(dat.lambda),toplot.AggG.(f_fig{i}),...
            'Linestyle','--','Linewidth',lw,'color',cmap(i,:)); hold on; % Phi during first month [line]
        plot(log10(dat.lambda),toplot.AggM.(f_fig{i})(:,1),...
            'Linewidth',lw,'color',cmap(i,:)); hold on; % Phi during last months [line]
        G(i) = scatter(log10(dat.lambda),toplot.AggG.(f_fig{i}),...
            'MarkerFaceColor','w','MarkerEdgeColor',cmap(i,:),...
            'SizeData',sz,'Linewidth',lw); hold on; % Phi during first month [scatter]
        M(i) = scatter(log10(dat.lambda),toplot.AggM.(f_fig{i}),...
            'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:),...
            'SizeData',sz,'Linewidth',lw); hold on; % Phi during last month [scatter]
        line([toplot.LG.(f_fig{i}) toplot.LG.(f_fig{i})],YL*YLIM(ifig,2),...
            'Linewidth',1.5*lw,'Linestyle',':','color',cmap(i,:)); % Diab. L during first month
        line([toplot.LM.(f_fig{i}) toplot.LM.(f_fig{i})],YL*YLIM(ifig,2),...
            'Linewidth',1.5*lw,'color',cmap(i,:)); % Diabatic. L during last month
    end; ylim(YLIM(ifig,:)); xlim(XLIM); grid on; % Enforce x and y bounds
    
    % Legend and Axes appearance
    L(ifig) = legend([G(FIE(1):FIE(end)) M(FIE(1):FIE(end))],LEG{2*(FIE(1)-1)+1:2*FIE(end)}); 
    set(L(ifig),'Location','Eastoutside','Fontsize',fz-2,'Interpreter','Latex');
    ax1=gca; set(ax1,'TickLabelInterpreter','Latex','Fontsize',fz); % Axes appearance
    ylabel(YLAB{ifig},'Fontsize',fz,'Interpreter','Latex'); % y label
    if ifig==3, G = gca; YTIK = G.YTickLabel; XTIK = G.XTickLabel;
        for ix = 1:numel(XTIK), XTIK{ix}=strcat('$10^{',XTIK{ix},'}$'); end
        for iy = 1:numel(YTIK), YTIK{iy}=strcat('$10^{',YTIK{iy},'}$'); end
        set(gca,'XTickLabel',XTIK); xlabel(XLAB,'Fontsize',fz,'Interpreter','Latex');
    end

end

%% 4. Adjust figure's position
ANOT = {'a)','b)','c)'};
for ifig = 1:3
    S(ifig).Position = [xleft+xleg (3-ifig)*yfig+(4-ifig)*ybot-(3-ifig)*ybot/3 xfig yfig];
    L(ifig).Position = [xfig+2*xleg (3-ifig)*yfig+(4-ifig)*ybot-(3-ifig)*ybot/3 1-xleft-xfig-xleg-dxleg yfig];
    A(ifig) = annotation('textbox',[S(ifig).Position(1) ...
        S(ifig).Position(2)+S(ifig).Position(4)-DelY DelX DelY],...
        'Fontsize',fz,'String',ANOT{ifig},'BackgroundColor','none',...
        'Interpreter','Latex','Linestyle','none','HorizontalAlignment','left'); hold on;
    axes('Position',S(ifig).Position,'Color','none',...
        'XAxisLocation','top','XLim',S(ifig).XLim,'YLim',S(ifig).YLim,...
        'XTickLabel','','YTickLabel','','Fontsize',fz,...
        'XTick',S(ifig).XTick,'YTick','','TickLabelInterpreter','Latex');
end

% Save Figure in PDF format
thisfile  = which(mfilename); basedir = thisfile(1:strfind(thisfile,mfilename)-1);
if strcmp(Lmet,'LAM')==1, name = 'A2'; else, name = 'A3'; end 
gcfsavepdf([basedir,'PDF_DATA',filesep,'Figure',name,'.pdf']); % Save plot
