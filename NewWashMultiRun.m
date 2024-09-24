warning('off')

%%

w0 = 2*pi*6.8944e-4; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
Q = 2e5; % Quality factor
kappa = I*w0^2; % Spring constance (N m/rad)
kb = 1.38064852e-23; % Boltzmann's constant (J/K)
T = 293; % Temperature (K)
thetaCalib = 3/300/8; % Autocollimator calibration (rad/(Diff/Sum))
m = 20e-3; % Mass (kg)
r = 2e-2; % Lever-arm (m)

Msun = 1.9891e30; % Mass of sun (kg)
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
Rsun = 149.6e9; % Radius from earth to sun (m)

aSun = G*Msun/Rsun^2; % Acceleration towards sun (m/s^2)
aEarth = 1.68e-2; % Acceleration towards center of Earth (m/s^2)
aGalaxy =  5e-11; % Acceleration towards dark matter at center of Galaxy (m/s^2)

TTFreq = 0.457120e-3; % Turn table frequency (Hz)

aDM = 1e-25/9.73e-19; % Torque to g_dm conversion for Be-Al

%% Injection Controls

injAmp = 20e-5*(r*m*aGalaxy);
inj1 = false; % Quadrature injection
inj2 = false; % Amplitude injection

%% Data loading

if (true)

    % Runs to load. Once per turntable cosine amplitude, sine amplitude, 
    % and misfit are calculated in NewWashAnalysis.m then loaded here
    runs = ["run6875Fits.mat" "run6891Fits.mat" "run6893Fits.mat" ...
        "run6895Fits.mat" "run6896Fits.mat" "run6897Fits.mat" "run6900Fits.mat" ...
        "run6903Fits.mat" "run6904Fits.mat"];

    timFitin =[];
    Cin = [];
    Sin = [];
    Uin = [];
    Pin = [];

    for f=1:length(runs)
        
        in = load(runs(f));

        % Chi-squared cut
        unCut = find(in.out(4,:).^2/(thermAmp^2)/sqrt(length(in.out(4,:))) < 5);

        timFitin = [timFitin in.out(1,unCut)];        
        Cin = [Cin detrend(in.out(2,unCut))];
        Sin = [Sin detrend(in.out(3,unCut))];
        Uin = [Uin in.out(4,unCut).^2/(thermAmp^2)/sqrt(length(in.out(4,:)))];
        % Pendulum flip after run6895
        if f<4
            Pin = [Pin in.out(4,unCut)*0+1];
        else 
            Pin = [Pin in.out(4,unCut)*0-1];
        end
    end
    
    % Moved time zero to midnight Jan. 1, 2024
    timFitin = mod(timFitin,31556926);

    %% 2-sigma cuts

    % Cosine cut
    unCut = find(and(Cin>prctile(Cin,5),Cin<prctile(Cin,95)));
    C = Cin(unCut);
    S = Sin(unCut);
    U = Uin(unCut);
    P = Pin(unCut);
    timFit = timFitin(unCut);

    % Sine cut
    unCut = find(and(S>prctile(S,5),S<prctile(S,95)));
    C = C(unCut);
    S = S(unCut);
    U = U(unCut);
    P = P(unCut);
    timFit = timFit(unCut);
    
    % Sampling frequency
    sampF = 1/(timFit(5)-timFit(4));

    % Loading in sun basis funtions outputted from sunVect.py
    rawSun=load('sunVectMin.out');
    sunSampF = 1/(rawSun(2,1)-rawSun(1,1))/3600/24;
    timSun=decimate(rawSun(:,1),floor(sunSampF/sampF));
    inSun=decimate(rawSun(:,2),floor(sunSampF/sampF));
    outSun=decimate(rawSun(:,3),floor(sunSampF/sampF));

    % Loading in galaxy basis funtions outputted from galVect.py
    rawGal=load('galVectMin.out');
    galSampF = 1/(rawGal(2,1)-rawGal(1,1))/3600/24;
    timGal=decimate(rawGal(:,1),floor(galSampF/sampF));
    inGal=decimate(rawGal(:,2),floor(galSampF/sampF));
    outGal=decimate(rawGal(:,3),floor(galSampF/sampF));
 
    % Quadrature injection
    if (inj1)
        for ind = 1:length(timFit)
            % Sync basis function and data
            galIndex = find(floor(timGal-timFit(ind)/24/3600)==0,1);
            C(ind) = C(ind) + sqrt(2)/2*P(ind)*injAmp*inGal(galIndex)+injAmp;
            S(ind) = S(ind) + sqrt(2)/2*P(ind)*injAmp*inGal(galIndex)+injAmp;
        end
    end
end

% Length of days
lenDays = ceil((timFit(end)-timFit(1))/24/3600);

% Calculate complex torque amplitude
torqFit = C+i*S;

% Fit parameters
fDay = 1/24/3600; % Daily frequency (Hz)
wDay = 2*pi*fDay; % Daily frequency (rad*Hz)
daySamples = floor(sampF/fDay); % Samples in a day
lenMin = 0.75*24*3600*TTFreq/2; % Minimum length of cut (samples)
numDaysFit = 2; % Length of cuts (days)

% Thermal noise
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq)); 

% Thermal noise circle
thermPhi = linspace(0,2*pi,100); 
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

%% Sun Fits

% Create vectors
cSun = [];
sSun = [];
timSunFit = [];
uSun = [];

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
    y = pol*detrend(abs(torqFit(indexCut)),'constant');

    % Amplitude injection
    if (inj2)
        for ind = 1:length(cut)
            % Sync basis function and data
            sunIndex = find(floor(timSun-cut(ind)/24/3600)==0,1);
            y(ind) = y(ind) + P(ind)*injAmp*inSun(sunIndex);
        end
    end

    if not(isempty(cut))   
        if (length(y)>=lenMin)

            % Sync basis function and data 
            sunIndex = [];
            for cutSun = cut
                sunIndex = [sunIndex find(floor(timSun-cutSun/24/3600)==0,1)];
            end

            % Design matrix
            x = [inSun(sunIndex) outSun(sunIndex)];

            if not(isempty(x))

                % Linear least squares fitting to basis functions
                a = inv(x'*x)*x'*y';

                % Append valid data points
                if (and(a(1)~=0,a(2)~=0))
                    cSun = [cSun a(1)];
                    sSun = [sSun a(2)];
                    uSun = [uSun std(a'*x'-y)];
                    timSunFit = [timSunFit mean(timFit(indexCut))];
                    longFit = [longFit a'*x'];
                    longTim = [longTim timFit(indexCut)];
                    longDat = [longDat y];
                    longU = [longU u];

                    % Add nans to long to allow gaps in plots
                    longFit = [longFit nan];
                    longTim = [longTim nan];
                    longDat = [longDat nan];
                    longU = [longU nan];
        
                end
            end
        end
    end
end

% Make complex sun torque 
torqSun = cSun+i*sSun;

% Mean and uncertainty of in-phase sun
ampSun = mean(cSun);
uncSun = std(cSun)/sqrt(length(cSun));

% Eotvos parameters
etaSun = ampSun/(r*m*aSun);
etaSunUnc = uncSun/(r*m*aSun);

% Display
disp(['Cosine Sun: ' num2str(mean(cSun)*1e15) ' fN m +- ' num2str(std(cSun)/sqrt(length(cSun))*1e15) ' fN m'])
disp(['Sine Sun: ' num2str(mean(sSun)*1e15) ' fN m +- ' num2str(std(sSun)/sqrt(length(sSun))*1e15) ' fN m'])
disp(['Amp Sun: ' num2str(ampSun*1e15) ' fN m +- ' num2str(uncSun*1e15) ' fN m'])
disp([' '])
disp(['Eta Sun: ' num2str(etaSun) ' +- ' num2str(etaSunUnc)])
disp([' '])

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
    y = pol*detrend(abs(torqFit(indexCut)),'constant');

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
                galIndex = [galIndex find(floor(timGal-cutGal/24/3600)==0,1)];
            end
            
            % Design matrix
            x = [inGal(galIndex) outGal(galIndex)];
            
            % Long basis functions for plotting 
            galLongIndex = find(and(timGal>=(index*numDaysFit+timFit(1)/24/3600),...
                timGal<=((index+1)*numDaysFit)+timFit(1)/24/3600));
            xLong = [inGal(galLongIndex) outGal(galLongIndex)];

            if not(isempty(x))

                % Linear least squares fitting to basis functions
                a = inv(x'*x)*x'*y';
            
                % Append valid data points
                if (and(a(1)~=0,a(2)~=0))
                    cGal = [cGal a(1)];
                    sGal = [sGal a(2)];
                    uGal = [uGal std(a'*x'-y)];
                    timGalFit = [timGalFit; timGal(galLongIndex)];
                    longFit = [longFit a'*xLong'];
                    longTim = [longTim timFit(indexCut)];
                    longDat = [longDat y];
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
uncGal = std(cGal)/sqrt(length(cGal));

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
figure(1)
l=plot(longTim/3600/24, longDat*1e15/(r*m), '.');
hold on
ll=plot(timGalFit,longFit*1e15/(r*m));
hold off
ylabel('Acceleration Amplitude (fm/s$^2$)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(ll,'LineWidth',2);
ylim([-30 30])
xlim([0.99*min(longTim/3600/24) 1.01*max(longTim/3600/24)])
legend('Data','Fit','Interpreter', 'latex')
grid on

% Sun fits
figure(2)
uncPhi = linspace(0,2*pi,100);
uncCirc = (std(cSun)/sqrt(length(cSun))*cos(thermPhi)+mean(cSun))...
        +i*(std(sSun)/sqrt(length(sSun))*sin(thermPhi)+mean(sSun));
l=plot(real(torqSun)*1e15/(r*m),imag(torqSun)*1e15/(r*m),'.');
hold on
ll=plot(mean(cSun)*1e15/(r*m),mean(sSun)*1e15/(r*m),'+',...
    real(uncCirc)*1e15/(r*m),imag(uncCirc)*1e15/(r*m), 'Color', [0.8500 0.3250 0.0980]);
hold off
ylabel('Out-of-Phase Sun (fm/s$^2$)','Interpreter', 'latex')
xlabel('In-Phase Sun (fm/s$^2$)','Interpreter', 'latex')
text(-2.5e-3,11,['$\eta_\odot$ = ' num2str(etaSun,2) ' $\pm$ ' num2str(etaSunUnc,2)],'Interpreter', 'latex','FontSize',16)
set(gca,'FontSize',16);
set(l,'MarkerSize',20);
set(ll,'LineWidth',2);
set(ll,'MarkerSize',12);
ylim([-15 15])
xlim([-15 15])
grid on

% Galactic DM fits
figure(3)
uncPhi = linspace(0,2*pi,100);
uncCirc = (std(cGal)/sqrt(length(cGal))*cos(thermPhi)+mean(cGal))...
        +i*(std(sGal)/sqrt(length(sGal))*sin(thermPhi)+mean(sGal));
l=plot(real(torqGal)*1e15/(r*m),imag(torqGal)*1e15/(r*m),'.');
hold on
ll=plot(mean(cGal)*1e15/(r*m),mean(sGal)*1e15/(r*m),'+',...
    real(uncCirc)*1e15/(r*m),imag(uncCirc)*1e15/(r*m), 'Color', [0.8500 0.3250 0.0980]);
hold off
ylabel('Out-of-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
xlabel('In-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
text(-2.5,12,['$\eta_{DM}$ = ' num2str(etaGalaxy,2) ' $\pm$ ' num2str(etaGalaxyUnc,2)],'Interpreter', 'latex','FontSize',16)
set(gca,'FontSize',16);
set(l,'MarkerSize',18);
set(ll,'LineWidth',2);
set(ll,'MarkerSize',12);
ylim([-15 15])
xlim([-15 15])
grid on

%% Save plots

if(false)
    fig2=figure(1);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_TimeSeries.pdf','-dpdf','-r1200')
    
    fig2=figure(2);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_SunFits.pdf','-dpdf','-r1200')

    fig2=figure(3);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_GalacticFits.pdf','-dpdf','-r1200')
end
