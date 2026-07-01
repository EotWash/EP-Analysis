warning('off')
clear
close all
%%

w0 = 2*pi*6.8567e-4; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
Q = 1.13e5; % Quality factor
kappa = I*w0^2; % Spring constance (N m/rad)
kb = 1.38064852e-23; % Boltzmann's constant (J/K)
T = 293; % Temperature (K)
thetaCalib = 0.0012523; % Autocollimator calibration (rad/(Diff/Sum))
m = 38.72e-3/2; % Mass (kg)
r = 3.77e-2/2; % Lever-arm (m)

Msun = 1.9891e30; % Mass of sun (kg)
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
Rsun = 149.6e9; % Radius from earth to sun (m)

aSun = G*Msun/Rsun^2; % Acceleration towards sun (m/s^2)
aEarth = 1.68e-2; % Acceleration towards center of Earth (m/s^2)
aGalaxy =  5e-11; % Acceleration towards dark matter at center of Galaxy (m/s^2)

sidDay = 86164.0905; % Sidereal day (s)
sidYear = 31558149.76; %Sideral year (s)

TTFreq = 0.4564e-3; % Turn table frequency (Hz)
filtPhase = 53.69*pi/180; % Low-pass filter phase correction (rad)
filtAt = 1/0.9958; % Low-pass filter amplitude correction (rad)

% Thermal noise
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq));

% Chi-squared threshold from distribution fit
thresh = 3.9794;

%% Data loading

if (true)

    % Runs to load. Once per turntable cosine amplitude, sine amplitude,
    % and misfit are calculated in NewWashAnalysis.m then loaded here
    runs = ["run6891Fits.mat" "run6893Fits.mat" ...
        "run6895Fits.mat" "run6896Fits.mat" "run6897Fits.mat" "run6900Fits.mat" ...
        "run6903Fits.mat" "run6904Fits.mat" "run6905Fits.mat" "run6923Fits.mat" ...
        "run6925Fits.mat" "run6926Fits.mat" "run6927Fits.mat" "run6930Fits.mat"...
        "run6931Fits.mat" "run6936Fits.mat" "run6939Fits.mat" "run6949Fits.mat"...
        "run6950Fits.mat" "run6954Fits.mat" "run6955Fits.mat" "run6956Fits.mat"...
        "run6958Fits.mat" "run6962Fits.mat" "run6964Fits.mat" "run6982Fits.mat"...
        "run6984Fits.mat", "run6985Fits.mat", "run6986Fits.mat", "run6987Fits.mat", "run6988Fits.mat"];
    runP = [1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1]; % 1 if 0 deg, 0 if 180 deg

    timFitin =[];
    Cin = [];
    Sin = [];
    Uin = [];
    Pin = [];
    Uraw = [];

    for f=1:length(runs)

        in = load("Fits\"+runs(f));

        % Chi-squared cut
        unCut = find(in.out(4,:)/thermAmp < thresh);

        % Moved time zero to midnight Jan. 1, 2024
        if f == 1
            timOffset = floor(in.out(1,1)/sidYear);
        end
        timFitin = [timFitin (in.out(1,unCut) - timOffset*sidYear)];

        Cin = [Cin detrend(in.out(2,unCut))];
        Sin = [Sin detrend(in.out(3,unCut))];
        Uin = [Uin in.out(4,unCut)/thermAmp];
        Uraw = [Uraw in.out(4,:)/thermAmp];

        % Pendulum flips
        if runP(f) %0 degrees
            Pin = [Pin in.out(4,unCut)*0+1];
        else %180 degrees
            Pin = [Pin in.out(4,unCut)*0-1];
        end
    end

    % 2-sigma cuts
    unCut = find(and(and(Cin>prctile(Cin,5),Cin<prctile(Cin,95)),...
        and(Sin>prctile(Sin,5),Sin<prctile(Sin,95))));
    C = Cin(unCut);
    S = Sin(unCut);
    U = Uin(unCut);
    P = Pin(unCut);
    timFit = mod(timFitin(unCut), sidYear);

    % Sampling frequency
    sampF = 1/(mode(round(diff(timFit),3)));

    % Loading in galaxy basis funtions outputted from galVect.py
    rawGal=load('Basis Functions\galVectMin.out');
    galSampF = 1/(rawGal(2,1)-rawGal(1,1))/sidDay;
    timGal=decimate(rawGal(:,1),floor(galSampF/sampF));
    inGal=detrend(decimate(rawGal(:,2),floor(galSampF/sampF)));
    outGal=detrend(decimate(rawGal(:,3),floor(galSampF/sampF)));

    % Daylight savings correction
    timGal = timGal - (timGal>307.042)/24 - (timGal<69.082)/24 - (timGal<433.041)/24 - (timGal>671.041)/24;

end

% Calculate complex torque amplitude
torqFit = filtAt*(cos(filtPhase)+i*sin(filtPhase))*P.*(C+i*S);

%% Best Fit

% Cut vectors
y = torqFit;

% Sync basis function and data
galIndex = [];
for cutGal = timFit
    [~,minI] = min(abs(timGal-cutGal/sidDay));
    galIndex = [galIndex minI];
end

% Design matrix
x = inGal(galIndex)+i*outGal(galIndex);

% Linear least squares fitting to basis functions
a = (x'*x)\x'*y';  

cGal = real(a(1));
sGal = imag(a(1));
uGal = std(a'*x'-y);
fitGal = a'*(inGal+i*outGal)';

% Complex galactic DM torque
torqGal = cGal+i*sGal;

%% Monte Carlo

mcN = 1e3; % Number of Monte Carlo samples

% Vector creation
cVect = [];
sVect = [];
Ntot = length(torqFit);
unVect = [];

for mcIndex = 1:mcN
    
    % Random point selection
    ranIndex = randi(length(torqFit),[length(torqFit),1]);

    % Cut vectors
    y = torqFit(ranIndex);
    
    % Sync basis function and data
    galIndex = [];
    for cutGal = timFit(ranIndex)
        [~,minI] = min(abs(timGal-cutGal/sidDay));
        galIndex = [galIndex minI];
    end
    
    % Design matrix
    x = [inGal(galIndex)+i*outGal(galIndex)];
    
    % Linear least squares fitting to basis functions
    a = (x'*x)\x'*y';  
    
    % Save values
    cVect = [cVect real(a(1))];
    sVect = [sVect imag(a(1))];
    unVect = [unVect std(cVect)];

end

% Statistics
ampGal = mean(cVect);
uncGal = std(cVect);

%%
% Eotvos parameters
etaGalaxy = ampGal/(r*m*aGalaxy);
etaGalaxyUnc = uncGal/(r*m*aGalaxy);

% Display
disp(['Cosine Galaxy: ' num2str(mean(cVect)/(r*m)*1e15) ' fm/s^2 +- ' num2str(std(cVect)/(r*m)/sqrt(length(cGal))*1e15) ' fm/s^2'])
disp(['Sine Galaxy: ' num2str(mean(sVect)/(r*m)*1e15) ' fm/s^2 +- ' num2str(std(sVect)/(r*m)/sqrt(length(sGal))*1e15) ' fm/s^2'])
disp(['Amp Galaxy: ' num2str(ampGal/(r*m)*1e15) ' fm/s^2 +- ' num2str(uncGal/(r*m)*1e15) ' fm/s^2'])
disp([' '])
disp(['Eta Galaxy: ' num2str(etaGalaxy) ' +- ' num2str(etaGalaxyUnc)])


%% Figures

% Time series
ylimit = 175;
textHeight = -85;

Rat = 9;
figure(1)
clf
set(gcf,'position',[100 100 1600 750]);
subplot(1,Rat,[1 Rat-1])
l=plot(mod(timFit,sidYear)/sidDay, torqFit*1e15/(r*m), '.');
hold on
ll=plot(mod(timGal,sidYear),fitGal*1e15/(r*m),[0 0],[0 0],...
    [40 40],[-ylimit ylimit],'k--', [114 114],[-ylimit ylimit],'k--', [183 183],[-ylimit ylimit],'k--', ...
    [211 211],[-ylimit ylimit],'k--', [284 284],[-ylimit ylimit],'k--');
text(10, textHeight, '$180^\circ$','Interpreter', 'latex','FontSize',16)
text(70, textHeight, '$0^\circ$','Interpreter', 'latex','FontSize',16)
text(140, textHeight, '$180^\circ$','Interpreter', 'latex','FontSize',16)
text(192, textHeight, '$0^\circ$','Interpreter', 'latex','FontSize',16)
text(243, textHeight, '$180^\circ$','Interpreter', 'latex','FontSize',16)
text(320, textHeight, '$0^\circ$','Interpreter', 'latex','FontSize',16)
patch([194.9 196.9 196.9 194.9],[-45 -45 45 45],'red','EdgeColor','k','FaceColor','none','LineWidth',1.5)
hold off
ylabel('Acceleration Amplitude (fm/s$^2$)','Interpreter', 'latex')
xlabel('Galactic Phase (sidereal days)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(ll,'LineWidth',1.5);
ylim([-ylimit ylimit])
xlim([0 365])
legend('Data', 'Fit','Fit $\times$ 100','Interpreter', 'latex')
grid on
%Hist subplot
subplot(1,Rat,Rat)
[n,xh] = hist(real(torqFit)*1e15/(r*m),15);
barh(xh,n,1);
ylim([-ylimit ylimit])
xlim([0 1.05*max(n)])
set(gca,'YTickLabel',[])
set(gca,'XGrid','off','YGrid','on')
set(gca, 'FontSize',16)

%Zoom Subplot
axes('position',[.325 .695 .15 .25])

zoomFirst = 184;
zoomEnd = 189;
zoomIndex = find(and(mod(timFit,sidYear)/sidDay>=zoomFirst,mod(timFit,sidYear)/sidDay<=zoomEnd));
zoomIndexGal = find(and(mod(timGal,sidYear)>=zoomFirst,mod(timGal,sidYear)<=zoomEnd));
l=plot(mod(timFit(zoomIndex),sidYear)/sidDay, real(torqFit(zoomIndex))*1e15/(r*m), '.');
hold on
ll=plot(mod(timGal(zoomIndexGal),sidYear),(real(fitGal(zoomIndexGal)))*1e15/(r*m),...
    mod(timGal(zoomIndexGal),sidYear),(real(fitGal(zoomIndexGal)))*1e15/(r*m)*1e2);
hold off
grid on
box on
set(gca, 'LineWidth',1.5)
xlim([zoomFirst zoomEnd])
ylim([-75 75])
set(l,'MarkerSize',16);
set(ll,'LineWidth',1.5);
set(gca, 'FontSize',14)

pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
width = pos(3);
height = pos(4);

annotation(gcf,'line',[0.4737 0.497],[0.6951 0.412],'LineWidth',1.7);
annotation(gcf,'line',[0.4747 0.496],[0.948 0.61],'LineWidth',1.7);

annotation(gcf,'rectangle',[0.135 0.135 0.675 0.1],'FaceColor',[0.90,0.90,0.90],'FaceAlpha',0);
annotation(gcf,'rectangle',[0.135 0.135 0.07 0.1],'FaceColor',[0.466666666666667 0.674509803921569 0.188235294117647],'FaceAlpha',0.15, 'LineStyle','none');
annotation(gcf,'rectangle',[0.205 0.135 0.14 0.1],'FaceColor',[0 0.282352941176471 0.470588235294118],'FaceAlpha',0.15, 'LineStyle','none');
annotation(gcf,'rectangle',[0.345 0.135 0.13 0.1],'FaceColor',[0.466666666666667 0.674509803921569 0.188235294117647],'FaceAlpha',0.15, 'LineStyle','none');
annotation(gcf,'rectangle',[0.475 0.135 0.051 0.1],'FaceColor',[0 0.282352941176471 0.470588235294118],'FaceAlpha',0.15, 'LineStyle','none');
annotation(gcf,'rectangle',[0.526 0.135 0.138 0.1],'FaceColor',[0.466666666666667 0.674509803921569 0.188235294117647],'FaceAlpha',0.15, 'LineStyle','none');
annotation(gcf,'rectangle',[0.664 0.135 0.146 0.1],'FaceColor',[0 0.282352941176471 0.470588235294118],'FaceAlpha',0.15, 'LineStyle','none');

annotation(gcf,'line',[0.15+10/width 0.15+10/width],[0.2 0.2-40/height],'LineWidth',3);
annotation(gcf,'ellipse',[0.15 0.2 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'ellipse',[0.15 0.2-40/height 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'arrow',[0.16 0.18],[0.2-10/height 0.2-10/height],'LineWidth',2,'Color',[0.301960784313725 0.745098039215686 0.933333333333333]);

annotation(gcf,'line',[0.26+10/width 0.26+10/width],[0.2 0.2-40/height],'LineWidth',3);
annotation(gcf,'ellipse',[0.26 0.2 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'ellipse',[0.26 0.2-40/height 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'arrow',[0.27 0.29],[0.2-10/height 0.2-10/height],'LineWidth',2,'Color',[0.301960784313725 0.745098039215686 0.933333333333333]);

annotation(gcf,'line',[0.4+10/width 0.4+10/width],[0.2 0.2-40/height],'LineWidth',3);
annotation(gcf,'ellipse',[0.4 0.2 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'ellipse',[0.4 0.2-40/height 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'arrow',[0.41 0.43],[0.2-10/height 0.2-10/height],'LineWidth',2,'Color',[0.301960784313725 0.745098039215686 0.933333333333333]);

annotation(gcf,'line',[0.49+10/width 0.49+10/width],[0.2 0.2-40/height],'LineWidth',3);
annotation(gcf,'ellipse',[0.49 0.2 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'ellipse',[0.49 0.2-40/height 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'arrow',[0.5 0.52],[0.2-10/height 0.2-10/height],'LineWidth',2,'Color',[0.301960784313725 0.745098039215686 0.933333333333333]);

annotation(gcf,'line',[0.59+10/width 0.59+10/width],[0.2 0.2-40/height],'LineWidth',3);
annotation(gcf,'ellipse',[0.59 0.2 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'ellipse',[0.59 0.2-40/height 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'arrow',[0.6 0.62],[0.2-10/height 0.2-10/height],'LineWidth',2,'Color',[0.301960784313725 0.745098039215686 0.933333333333333]);

annotation(gcf,'line',[0.73+10/width 0.73+10/width],[0.2 0.2-40/height],'LineWidth',3);
annotation(gcf,'ellipse',[0.73 0.2 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'ellipse',[0.73 0.2-40/height 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'arrow',[0.74 0.76],[0.2-10/height 0.2-10/height],'LineWidth',2,'Color',[0.301960784313725 0.745098039215686 0.933333333333333]);
%%
% Galactic DM fits
figure(3)
clf
set(gcf,'position',[100 100 900 900]);
tiledlayout(4,4)
nexttile([1,3])
x = [-2.5 -2.0 -1.5 -1.0 -.5 0 .5 1.0 1.5 2.0 2.5];

labelH = 320;

[n,x] = hist(cVect*1e15/(r*m), x);
bar(x,n,1);
hold on
text(-2.4,labelH,['$\mu_{in}$ = ' num2str(mean(cVect)/(r*m)*1e15,1) ' fm/s$^2$'],'Interpreter', 'latex','FontSize',14)
text(1.1,labelH,['$\sigma_{in}$ = ' num2str(std(cVect)/(r*m)*1e15,2) ' fm/s$^2$'],'Interpreter', 'latex','FontSize',14)
hold off
xlim([min(x) max(x)])
ylim([0 max(n)*1.2])
set(gca,'XTickLabel',[])
set(gca,'XTick',x)
xticks(x)
set(gca,'XGrid','on','YGrid','off')
set(gca,'FontSize',16);
nexttile(5,[3,3])
x2 = linspace(-5, 5, 20);
[N,xe,ye]=histcounts2(cVect*1e15/(r*m),sVect*1e15/(r*m),x2+0.25,x2+0.25);
c=contourf(xe(1:end-1)+0.25,ye(1:end-1)+0.25,N',linspace(1, max(max(N)),4),'edgecolor','none');
view(0,90)
colormap(sky);
shading(gca,'interp')
ylabel('Out-of-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
xlabel('In-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
set(gca,'FontSize',16);
ylim([min(x) max(x)])
xlim([min(x) max(x)])
set(gca,'XTick',x)
set(gca,'YTick',x)
xticks(x)
yticks(x)
grid on

nexttile(8,[3,1])
[n,x] = hist(sVect*1e15/(r*m),x);
barh(x,n,1);
hold on
text(labelH,2.4,['$\mu_{out}$ = ' num2str(mean(sVect)/(r*m)*1e15,1) ' fm/s$^2$'],'Interpreter', 'latex','FontSize',14,'Rotation',-90)
text(labelH,-1.1,['$\sigma_{out}$ = ' num2str(std(sVect)/(r*m)*1e15,1) ' fm/s$^2$'],'Interpreter', 'latex','FontSize',14,'Rotation',-90)
hold off
ylim([min(x) max(x)])
xlim([0 max(n)*1.2])
set(gca,'YTickLabel',[])
set(gca,'XGrid','off','YGrid','on')
set(gca,'YTick',x)
yticks(x)
set(gca,'FontSize',16);


%% Save plots

if(true)
    fig2=figure(1);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_TimeSeries.pdf','-dpdf','-r1200')

    fig2=figure(3);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_Histogram.pdf','-dpdf','-r1200')
end
