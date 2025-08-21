warning('off')

%%

w0 = 2*pi*6.8567e-4; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
% Q = 2.89e5; % Quality factor
Q = 1.13e5;
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
% aGalaxy =  9.7e-11; % Acceleration towards dark matter at center of Galaxy (m/s^2)

sidDay = 86164.0905; % Sidereal day (s)

TTFreq = 0.457120e-3; % Turn table frequency (Hz)

% Thermal noise
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq)); 

% %% Injection Controls
% inj = [];
% injEx = [];
% injUnc = [];

% inj = [inj injAmp/(r*m*aGalaxy)];
% injEx = [injEx etaGalaxy];
% injUnc = [injUnc etaGalaxyUnc];

injAmp = 10e-5*(r*m*aGalaxy);
inj1 = false; % Quadrature injection
inj2 = false; % Amplitude injection

% Chi-squared threshold
thresh = 7;

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
        "run6958Fits.mat" "run6962Fits.mat" "run6964Fits.mat"];

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
        if f<10
            % 2024 Runs
            timFitin = [timFitin mod(in.out(1,unCut),31556926)];
        else
            % 2025 Runs
            timFitin = [timFitin mod(in.out(1,unCut),31556926)+31556926];
        end
        Cin = [Cin detrend(in.out(2,unCut))];
        Sin = [Sin detrend(in.out(3,unCut))];
        Uin = [Uin in.out(4,unCut)/thermAmp];
        Uraw = [Uraw in.out(4,:)/thermAmp];
        
        % Pendulum flips
        if or(f<3,and(f>13,f<20))
            Pin = [Pin in.out(4,unCut)*0+1];
        else 
            Pin = [Pin in.out(4,unCut)*0-1];
        end
    end

    %% 2-sigma cuts

    % Cosine cut
    unCut = find(and(and(Cin>prctile(Cin,5),Cin<prctile(Cin,95)),...
        and(Sin>prctile(Sin,5),Sin<prctile(Sin,95))));
    C = Cin(unCut);
    S = Sin(unCut);
    U = Uin(unCut);
    P = Pin(unCut);
    timFit = timFitin(unCut);
    
    % Sampling frequency
    sampF = 1/(timFit(5)-timFit(4));

    % Loading in galaxy basis funtions outputted from galVect.py
    rawGal=load('Basis Functions\galVectMin.out');
    galSampF = 1/(rawGal(2,1)-rawGal(1,1))/3600/24;
    timGal=decimate(rawGal(:,1),floor(galSampF/sampF));
    inGal=(decimate(rawGal(:,2),floor(galSampF/sampF)));
    outGal=detrend(decimate(rawGal(:,3),floor(galSampF/sampF)));

    timGal = timGal - (timGal>307.042)/24 - (timGal<69.082)/24 - (timGal<433.041)/24 - (timGal>671.041)/24;

end

% Length of days
lenDays = ceil((timFit(end)-timFit(1))/24/3600);

% Calculate complex torque amplitude
torqFit = P.*(C+i*S);

% Quadrature injection
if (inj1)    
    torqFit = injAmp*inGal'+abs(injAmp);
    timFit = timGal'*24*3600;
    P = ones(length(inGal),1);
    U = 0*P';
end


% Fit parameters
fDay = 1/sidDay; % Daily frequency (Hz)
wDay = 2*pi*fDay; % Daily frequency (rad*Hz)
daySamples = floor(sampF/fDay); % Samples in a day
numDaysFit = 2; % Length of cuts (days)
lenMin = 0*numDaysFit; % Minimum length of cut (samples)

% Thermal noise circle
thermPhi = linspace(0,2*pi,100); 
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

% Remove insanely large chi-squareds
Uraw = Uraw(find(Uraw<20));

%% Galaxy Fits

% Create vectors
cGal = [];
sGal = [];
timGalFit = [];
uGal = [];

% Create long vectors for plotting
longFit = [];
longTim = [];
longDat = [];
longU = [];

for index = 0:lenDays/numDaysFit

    % Find cut indices
    indexCut = find(floor((timFit-timFit(1))/24/3600/numDaysFit)==index);

    % Cut vectors
    cut = timFit(indexCut);
    pol = sign(mean(P(indexCut)));
    u = U(indexCut);
    y = torqFit(indexCut)-mean(torqFit(indexCut));

    % Amplitude injection
    if (inj2)
        for ind = 1:length(cut)
            % Sync basis function and data
            galIndex = find(floor(timGal-cut(ind)/24/3600)==0,1);
            y(ind) = y(ind) + P(ind)*injAmp*inGal(galIndex);
        end
    end

    if not(isempty(cut))   
        if (length(y)>=lenMin)

            % Sync basis function and data 
            galIndex = [];
            for cutGal = cut
%                 galIndex = [galIndex find(floor(timGal-cutGal/24/3600)==0,1)];
                [m,minI] = min(abs(timGal-cutGal/24/3600));
                galIndex = [galIndex minI];
            end
            
            % Design matrix
            x = [inGal(galIndex)+i*outGal(galIndex)];
            
            % Long basis functions for plotting 
            galLongIndex = find(and(timGal>=(index*numDaysFit+timFit(1)/24/3600),...
                timGal<=((index+1)*numDaysFit)+timFit(1)/24/3600));
            xLong = [inGal(galLongIndex)+i*outGal(galLongIndex)];

            if not(isempty(x))

                % Linear least squares fitting to basis functions
                a = inv(x'*x)*x'*y';
            
                % Append valid data points
                if ((a(1)~=0))
                    cGal = [cGal real(a(1))];
                    sGal = [sGal imag(a(1))];
                    uGal = [uGal std(a'*x'-y)];
                    timGalFit = [timGalFit; timGal(galLongIndex)];
                    longFit = [longFit (a'*xLong')-mean(a'*xLong')];
                    longTim = [longTim timFit(indexCut)];
                    longDat = [longDat (real(y)+imag(y))/2];
                    longU = [longU u];
        
                    % Add nans to long to allow gaps in plots
                    longFit = [longFit nan];
                    longTim = [longTim nan];
                    longDat = [longDat nan];
                    longU = [longU nan];
                    timGalFit = [timGalFit; nan];
        
                end
            end
        end
    end
end

% Complex galactic DM torque
torqGal = cGal+i*sGal;

% Mean and uncertainty of in-phase torque
ampGal = mean(cGal);
uncGal = std(cGal)/sqrt(length(sGal));

%% Uncertainty vs time calc
cGalA = (cGal);
uncT = 1:length(cGal);
uncGalT = [];
for index=uncT
    cGalT = cGalA(1:index);
    uncGalT = [uncGalT; std(cGalT)/sqrt(length(cGalT))/(r*m*aGalaxy)];

end

%%

% Eotvos parameters
etaGalaxy = ampGal/(r*m*aGalaxy);
etaGalaxyUnc = uncGal/(r*m*aGalaxy);

% Display
disp(['Cosine Galaxy: ' num2str(mean(cGal)/(r*m)*1e15) ' fm/s^2 +- ' num2str(std(cGal)/(r*m)/sqrt(length(cGal))*1e15) ' fm/s^2'])
disp(['Sine Galaxy: ' num2str(mean(sGal)/(r*m)*1e15) ' fm/s^2 +- ' num2str(std(sGal)/(r*m)/sqrt(length(sGal))*1e15) ' fm/s^2'])
disp(['Amp Galaxy: ' num2str(ampGal/(r*m)*1e15) ' fm/s^2 +- ' num2str(uncGal/(r*m)*1e15) ' fm/s^2'])
disp([' '])
disp(['Eta Galaxy: ' num2str(etaGalaxy) ' +- ' num2str(etaGalaxyUnc)])

%% Figures


% Time series
Rat = 9;
figure(1)
set(gcf,'position',[100 100 1600 750]);
subplot(1,Rat,[1 Rat-1])
l=plot(longTim/3600/24, longDat*1e15/(r*m), '.');
hold on
ll=plot(timGalFit,longFit*1e15/(r*m), [213 213],[-80 80],'k--', [420 420],[-80 80],'k--', [486 486],[-80 80],'k--');
% patch([275 382 382 275], [-80 -80 80 80], [.5 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5)
% text(300, 2, 'Hardware Failures','Interpreter', 'latex','FontSize',16)
% text(200, -70, '0$^\circ$','Interpreter', 'latex','FontSize',16)
% text(235, -70, '180$^\circ$','Interpreter', 'latex','FontSize',16)
% text(393, -70, '180$^\circ$','Interpreter', 'latex','FontSize',16)
% text(455, -70, '0$^\circ$','Interpreter', 'latex','FontSize',16)
% text(492, -70, '180$^\circ$','Interpreter', 'latex','FontSize',16)
hold off
ylabel('Acceleration Amplitude (fm/s$^2$)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(ll,'LineWidth',2);
ylim([-80 80])
xlim([0.96*min(longTim/3600/24) 1.01*max(longTim/3600/24)])
legend('Data','Fit','Interpreter', 'latex')
grid on

subplot(1,Rat,Rat)
[n,x] = hist(longDat*1e15/(r*m),15);
barh(x,n,1);
% hold on
% xTherm = linspace(-1.5*max(x),1.5*max(x));
% plot(max(n)*exp(-(xTherm.^2)/((thermAmp/sqrt(2)*1e15/r/m)^2)),xTherm,'LineWidth',2)
% hold off
ylim([-80 80])
xlim([0 1.01*max(n)])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'XGrid','off','YGrid','on')
% legend('Data','Thermal Noise','Interpreter', 'latex')

axes('position',[.32 .6 .15 .25])
zoomIndex = find(and(longTim/3600/24>=384.8,longTim/3600/24<=386.8));
l=plot(longTim(zoomIndex)/3600/24, longDat(zoomIndex)*1e15/(r*m), '.');
hold on
zoomIndex2 = find(and(timGalFit>=384.8,timGalFit<=386.8));
ll=plot(timGalFit(zoomIndex2),longFit(zoomIndex2)*1e15/(r*m));
hold off
grid on
box on
set(gca, 'LineWidth',1.5)
xlim([384.8 386.8])
ylim([-40 40])
set(l,'MarkerSize',8);
set(ll,'LineWidth',1.5);

annotation(gcf,'line',[0.47 0.5],[0.6 0.31],'LineWidth',1.5);
annotation(gcf,'line',[0.47 0.5],[0.85 0.72],'LineWidth',1.5);
annotation(gcf,'rectangle',[0.5 0.31 0.005 0.41],'LineWidth',1.5);

pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
width = pos(3);
height = pos(4);

annotation(gcf,'ellipse',[0.14 0.2 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'ellipse',[0.14 0.2-40/height 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'ellipse',[0.14+40/width 0.2 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'ellipse',[0.14+40/width 0.2-40/height 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);

annotation(gcf,'ellipse',[0.35 0.2 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'ellipse',[0.35 0.2-40/height 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'ellipse',[0.35+40/width 0.2 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'ellipse',[0.35+40/width 0.2-40/height 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);

annotation(gcf,'ellipse',[0.6 0.2 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'ellipse',[0.6 0.2-40/height 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'ellipse',[0.6+40/width 0.2 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'ellipse',[0.6+40/width 0.2-40/height 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);

annotation(gcf,'ellipse',[0.72 0.2 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'ellipse',[0.72 0.2-40/height 20/width 20/height],'FaceColor',[0.85,0.33,0.10]);
annotation(gcf,'ellipse',[0.72+40/width 0.2 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);
annotation(gcf,'ellipse',[0.72+40/width 0.2-40/height 20/width 20/height],'FaceColor',[0.93,0.69,0.13]);


%%
% Galactic DM fits
figure(3)
uncPhi = linspace(0,2*pi,100);
uncCirc = (2*std(cGal)/sqrt(length(cGal))*cos(thermPhi)+mean(cGal))...
        +i*(2*std(sGal)/sqrt(length(sGal))*sin(thermPhi)+mean(sGal));
tiledlayout(4,4)
nexttile([1,3])
x = [-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30];
[n,x] = hist(real(torqGal)*1e15/(r*m), x);
bar(x,n,1);
xlim([-30 30])
% xlim([0 150])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'XTick',x)
set(gca,'XGrid','on','YGrid','off')

nexttile(5,[3,3])
l=plot(real(torqGal)*1e15/(r*m),imag(torqGal)*1e15/(r*m),'.');
hold on
ll=plot(mean(cGal)*1e15/(r*m),mean(sGal)*1e15/(r*m),'+',...
    real(uncCirc)*1e15/(r*m),imag(uncCirc)*1e15/(r*m), 'Color', [0.8500 0.3250 0.0980]);
hold off
ylabel('Out-of-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
xlabel('In-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
text(-10,17,['$\eta_{DM}$ = (' num2str(1e5*etaGalaxy,4) ' $\pm$ ' num2str(1e5*etaGalaxyUnc,4) ') $\times\ 10^{-5}$'],'Interpreter', 'latex','FontSize',16)
set(gca,'FontSize',16);
set(l,'MarkerSize',18);
set(ll,'LineWidth',2);
set(ll,'MarkerSize',12);
ylim([-30 30])
xlim([-30 30])
set(gca,'XTick',x)
set(gca,'YTick',x)

grid on

nexttile(8,[3,1])
[n,x] = hist(imag(torqGal)*1e15/(r*m),x);
barh(x,n,1);
ylim([-30 30])
% xlim([0 150])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'XGrid','off','YGrid','on')
set(gca,'YTick',x)

figure(5)
[NR,XR] = hist(Uraw,100);
[NU,XU] = hist(U,XR);
bar(XR,NR,'FaceAlpha',0.5)
hold on
bar(XU,NU,'FaceAlpha',0.5)
hold off
% xlim([-4 4]*1e-25)
xlabel('$\chi^2$','Interpreter', 'latex')
ylabel('Number','Interpreter', 'latex')
legend('Without Data Quality Cut','With Data Quality Cut','Interpreter', 'latex')
set(gca,'FontSize',16);
grid on
%%
figure(6)
ll = plot(uncT*2, uncGalT, linspace(0,200), thermAmp/(r*m*aGalaxy)/sqrt(24*3600*TTFreq/numDaysFit/3)./sqrt(linspace(0,200)),...
    [0 200], [1.3e-5 1.3e-5],'--');
% xlim([-4 4]*1e-25)
xlabel('Number of Days','Interpreter', 'latex')
ylabel('Uncertainty','Interpreter', 'latex')
legend('Data', 'Thermal Noise','$\eta < 1.3 \times 10^{-5}$','Interpreter', 'latex')
set(gca,'FontSize',16);
set(ll,'LineWidth',2);
ylim([0 4e-5])
grid on

%% Save plots

if(false)
    fig2=figure(1);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_TimeSeries.pdf','-dpdf','-r1200')
    
    fig2=figure(3);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_GalacticFits.pdf','-dpdf','-r1200')

    fig2=figure(6);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_UncVsTime.pdf','-dpdf','-r1200')
end
