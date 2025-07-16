warning('off')

%%

w0 = 2*pi*6.8567e-4; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
Q = 2.89e5; % Quality factor
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

TTFreq = 0.457120e-3; % Turn table frequency (Hz)

% Thermal noise
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq)); 

%% Injection Controls

injAmp = 200e-5*(r*m*aSun);
inj1 = false; % Quadrature injection
inj2 = false; % Amplitude injection

%% Data loading

if (true)

    % Runs to load. Once per turntable cosine amplitude, sine amplitude, 
    % and misfit are calculated in NewWashAnalysis.m then loaded here
    runs = ["run6891Fits.mat" "run6893Fits.mat" ...
       "run6895Fits.mat" "run6896Fits.mat" "run6897Fits.mat" "run6900Fits.mat" ...
        "run6903Fits.mat" "run6904Fits.mat" "run6905Fits.mat" "run6923Fits.mat" ...
        "run6925Fits.mat" "run6926Fits.mat" "run6927Fits.mat" "run6930Fits.mat"...
        "run6931Fits.mat" "run6936Fits.mat" "run6939Fits.mat"];

    timFitin =[];
    Cin = [];
    Sin = [];
    Uin = [];
    Pin = [];
    Uraw = [];

    for f=1:length(runs)
        
        in = load(runs(f));        

        % Chi-squared cut        
        unCut = find(in.out(4,:).^2/(thermAmp^2)/sqrt(length(in.out(4,:))) < 4);
        
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
        Uin = [Uin in.out(4,unCut).^2/(thermAmp^2)/sqrt(length(in.out(4,:)))];
        Uraw = [Uraw in.out(4,:).^2/(thermAmp^2)/sqrt(length(in.out(4,:)))];
        
        % Pendulum flip after run6895
        if or(f<3,f>13)
            Pin = [Pin in.out(4,unCut)*0+1];
        else 
            Pin = [Pin in.out(4,unCut)*0-1];
        end
        if f==2
            timFlip = mod(in.out(1,end),31556926);
        end
    end

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

    timSun = timSun - (timSun>307.042)/24 - (timSun<69.082)/24 - (timSun<433.041)/24 - (timSun>671.041)/24;

    timSun = [timSun; timSun+365.25];
    inSun = [inSun; inSun];
    outSun = [outSun; outSun];
 
    % Quadrature injection
    if (inj1)
        for ind = 1:length(timFit)
            % Sync basis function and data
            SunIndex = find(floor(timSun-timFit(ind)/24/3600)==0,1);
            C(ind) = C(ind) + sqrt(2)/2*P(ind)*injAmp*inSun(SunIndex)+injAmp;
            S(ind) = S(ind) + sqrt(2)/2*P(ind)*injAmp*inSun(SunIndex)+injAmp;
        end
    end
end

% Length of days
lenDays = ceil((timFit(end)-timFit(1))/24/3600);

% Calculate complex torque amplitude
torqFit = P.*(C+i*S);

% Fit parameters
fDay = 1/24/3600; % Daily frequency (Hz)
wDay = 2*pi*fDay; % Daily frequency (rad*Hz)
daySamples = floor(sampF/fDay); % Samples in a day
numDaysFit = 2; % Length of cuts (days)
lenMin = 5*numDaysFit; % Minimum length of cut (samples)

% Thermal noise circle
thermPhi = linspace(0,2*pi,100); 
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

% Remove insanely large chi-squareds
Uraw = Uraw(find(Uraw<20));

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
    y = torqFit(indexCut)-mean(torqFit(indexCut));

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
            x = [inSun(sunIndex)+i*outSun(sunIndex)];

            % Long basis functions for plotting 
            sunLongIndex = find(and(timSun>=(index*numDaysFit+timFit(1)/24/3600),...
                timSun<=((index+1)*numDaysFit)+timFit(1)/24/3600));
            xLong = [inSun(sunLongIndex)+i*outSun(sunLongIndex)];

            if not(isempty(x))

                % Linear least squares fitting to basis functions
                a = inv(x'*x)*x'*y';

                % Append valid data points
                if (a(1)~=0)
                    cSun = [cSun real(a(1))];
                    sSun = [sSun imag(a(1))];
                    uSun = [uSun std(a'*x'-y)];
                    timSunFit = [timSunFit; timSun(sunLongIndex)];
                    longFit = [longFit a'*xLong'];
                    longTim = [longTim timFit(indexCut)];
                    longDat = [longDat (real(y)+imag(y))/2];
                    longU = [longU u];

                    % Add nans to long to allow gaps in plots
                    longFit = [longFit nan];
                    longTim = [longTim nan];
                    longDat = [longDat nan];
                    longU = [longU nan];
                    timSunFit = [timSunFit; nan];
        
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
disp(['Cosine Sun: ' num2str(mean(cSun/r/m)*1e15) ' fm/s^2+- ' num2str(std(cSun/r/m)/sqrt(length(cSun))*1e15) ' fm/s^2'])
disp(['Sine Sun: ' num2str(mean(sSun/r/m)*1e15) ' fm/s^2 +- ' num2str(std(sSun/r/m)/sqrt(length(sSun))*1e15) ' fm/s^2'])
disp(['Amp Sun: ' num2str(ampSun/r/m*1e15) ' fm/s^2 +- ' num2str(uncSun/r/m*1e15) ' fm/s^2'])
disp([' '])
disp(['Eta Sun: ' num2str(etaSun) ' +- ' num2str(etaSunUnc)])
disp([' '])

%% Figures

% Time series
Rat = 9;
figure(1)
subplot(1,Rat,[1 Rat-1])
l=plot(longTim/3600/24, longDat*1e15/(r*m), '.');
hold on
ll=plot(timSunFit,longFit*1e15/(r*m), [213 213],[-45 45],'k--', [420 420],[-45 45],'k--');
patch([275 384 384 275], [-45 -45 45 45], [.5 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5)
text(315, 2, 'Hardware Failures','Interpreter', 'latex','FontSize',16)
text(200, 32, '0$^\circ$','Interpreter', 'latex','FontSize',16)
text(235, 32, '180$^\circ$','Interpreter', 'latex','FontSize',16)
text(392, 32, '180$^\circ$','Interpreter', 'latex','FontSize',16)
text(435, 32, '0$^\circ$','Interpreter', 'latex','FontSize',16)
hold off
ylabel('Acceleration Amplitude (fm/s$^2$)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(ll,'LineWidth',2);
ylim([-45 45])
xlim([0.99*min(longTim/3600/24) 1.01*max(longTim/3600/24)])
legend('Data','Fit','Interpreter', 'latex')
grid on

subplot(1,Rat,Rat)
[n,x] = hist(longDat*1e15/(r*m));
barh(x,n,1);
hold on
xTherm = linspace(-1.5*max(x),1.5*max(x));
plot(max(n)*exp(-(xTherm.^2)/(thermAmp*1e15/r/m/sqrt(2))^2),xTherm, '--','LineWidth',4)
hold off
ylim([-30 30])
xlim([0 1.01*max(n)])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'XGrid','off','YGrid','on')
% legend('Data','Thermal Noise','Interpreter', 'latex')

% Sun fits
figure(3)
uncPhi = linspace(0,2*pi,100);
uncCirc = (std(cSun)/sqrt(length(cSun))*cos(thermPhi)+mean(cSun))...
        +i*(std(sSun)/sqrt(length(sSun))*sin(thermPhi)+mean(sSun));
tiledlayout(4,4)
nexttile([1,3])
x = linspace(-15,15,10);
[n,x] = hist(real(torqSun)*1e15/(r*m), x);
bar(x,n,1);
xlim([-15 15])
% xlim([0 150])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'XTick',[-10 -5 0 5 10])
set(gca,'XGrid','on','YGrid','off')

nexttile(5,[3,3])
l=plot(real(torqSun)*1e15/(r*m),imag(torqSun)*1e15/(r*m),'.');
hold on
ll=plot(mean(cSun)*1e15/(r*m),mean(sSun)*1e15/(r*m),'+',...
    real(uncCirc)*1e15/(r*m),imag(uncCirc)*1e15/(r*m), 'Color', [0.8500 0.3250 0.0980]);
hold off
ylabel('Out-of-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
xlabel('In-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
% text(-4,9,['$\eta_\odot$ = ' num2str(etaSun,2) ' $\pm$ ' num2str(etaSunUnc,2)],'Interpreter', 'latex','FontSize',16)
set(gca,'FontSize',16);
set(l,'MarkerSize',18);
set(ll,'LineWidth',2);
set(ll,'MarkerSize',12);
ylim([-15 15])
xlim([-15 15])
grid on

nexttile(8,[3,1])
[n,x] = hist(imag(torqSun)*1e15/(r*m),x);
barh(x,n,1);
ylim([-15 15])
% xlim([0 150])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'XGrid','off','YGrid','on')
set(gca,'YTick',[-10 -5 0 5 10])

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
%% Save plots

if(false)
    fig2=figure(1);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_SunTimeSeries.pdf','-dpdf','-r1200')
    
    fig2=figure(3);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_SunFits.pdf','-dpdf','-r1200')
    
end
