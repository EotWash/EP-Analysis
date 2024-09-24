warning('off')

%%

w0 = 2*pi*6.8944e-4;
I = 3.78e-5;
Q = 2e5;
kappa = I*w0^2;
% kappa = 7.129e-10; % Theoretical 
kb = 1.38064852e-23;
T = 293;
thetaCalib = 3/300/8;
m = 20e-3;
r = 2e-2;

Msun = 1.9891e30; % Mass of sun (kg)
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
Rsun = 149.6e9; % Radius from earth to sun (m)

aSun = G*Msun/Rsun^2; % Acceleration towards sun (m/s^2)
aEarth = 1.68e-2;
aGalaxy =  5e-11;

% TTFreq = mean(diff(TTAngle/360))*sampF;
TTFreq = 0.457120e-3;

aDM = 1e-25/9.73e-19; % Torque to g_dm conversion for Be-Al

%% Thermal Circle Calculations

thermPhi = linspace(0,2*pi,100);
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq));
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

%%

injAmp = 20e-5*(r*m*aGalaxy);
inj1 = false; % Quadrature injection
inj2 = false; % Amplitude injection

%% Data loading

if (true)

    runs = ["run6875Fits.mat" "run6891Fits.mat" "run6893Fits.mat" ...
        "run6895Fits.mat" "run6896Fits.mat" "run6897Fits.mat" "run6900Fits.mat" ...
        "run6903Fits.mat" "run6904Fits.mat"];
%     runs = ["run6891Fits.mat"];
    timFitin =[];
    Cin = [];
    Sin = [];
    Uin = [];
    Pin = [];
    for f=1:length(runs)
        in = load(runs(f));

        unCut = find(in.out(4,:).^2/(thermAmp^2)/sqrt(length(in.out(4,:))) < 5);

        timFitin = [timFitin in.out(1,unCut)];
        Cin = [Cin detrend(in.out(2,unCut))];
        Sin = [Sin detrend(in.out(3,unCut))];
        Uin = [Uin in.out(4,unCut).^2/(thermAmp^2)/sqrt(length(in.out(4,:)))];
        if f<4
            Pin = [Pin in.out(4,unCut)*0+1];
        else 
            Pin = [Pin in.out(4,unCut)*0-1];
        end
    end
    
    timFitin = mod(timFitin,31556926);

    inDEP = load('DarkEPData.csv');
    fDEP = inDEP(:,1);
    aDEP = inDEP(:,2);
   

    %% Uncertainty cuts

    sig = std(Cin);

    unCut = find(and(Cin>prctile(Cin,5),Cin<prctile(Cin,95)));
    C = Cin(unCut);
    S = Sin(unCut);
    U = Uin(unCut);
    P = Pin(unCut);
    timFit = timFitin(unCut);

    unCut = find(and(S>prctile(S,5),S<prctile(S,95)));
    C = C(unCut);
    S = S(unCut);
    U = U(unCut);
    P = P(unCut);
    timFit = timFit(unCut);
    
    sampF = 1/(timFit(5)-timFit(4));

    % Loading in basis funtions outputted from sunVect.py
    rawSun=load('sunVectMin.out');
    sunSampF = 1/(rawSun(2,1)-rawSun(1,1))/3600/24;
    timSun=decimate(rawSun(:,1),floor(sunSampF/sampF));
    inSun=decimate(rawSun(:,2),floor(sunSampF/sampF));
    outSun=decimate(rawSun(:,3),floor(sunSampF/sampF));

    % Loading in basis funtions outputted from sunVect.py
    rawGal=load('galVectMin.out');
    galSampF = 1/(rawGal(2,1)-rawGal(1,1))/3600/24;
    timGal=decimate(rawGal(:,1),floor(galSampF/sampF));
    inGal=decimate(rawGal(:,2),floor(galSampF/sampF));
    outGal=decimate(rawGal(:,3),floor(galSampF/sampF));
 
    if (inj1)
        for ind = 1:length(timFit)
            % Sync basis function and data
            galIndex = find(floor(timGal-timFit(ind)/24/3600)==0,1);
            C(ind) = C(ind) + sqrt(2)/2*P(ind)*injAmp*inGal(galIndex)+injAmp;
            S(ind) = S(ind) + sqrt(2)/2*P(ind)*injAmp*inGal(galIndex)+injAmp;
        end
    end
end

lenDays = ceil((timFit(end)-timFit(1))/24/3600);
%%

torqFit = C+i*S;

amp = sqrt(mean(C)^2+mean(S)^2);
unc = 1/amp*sqrt(mean(C)^2*std(C)^2+mean(S)^2*std(S)^2)/sqrt(length(C));

disp(['Cosine Torque: ' num2str(mean(C)*1e15) ' fN m +- ' num2str(std(C)/sqrt(length(C))*1e15) ' fN m'])
disp(['Sine Torque: ' num2str(mean(S)*1e15) ' fN m +- ' num2str(std(S)/sqrt(length(S))*1e15) ' fN m'])
disp(['Amp Torque: ' num2str(amp*1e15) ' fN m +- ' num2str(unc*1e15) ' fN m'])
disp([' '])



%% Sun Fits

cSun = [];
sSun = [];
timSunFit = [];
uSun = [];

longFit = [];
longTim = [];
longDat = [];
longU = [];

fDay = 1/24/3600;
wDay = 2*pi*fDay;
daySamples = floor(sampF/fDay);
lenMin = 0.75*24*3600*TTFreq/2;
numDaysFit = 2;

for index = 0:lenDays/numDaysFit

    indexCut = find(floor((timFit-timFit(1))/24/3600/numDaysFit)==index);
    cut = timFit(indexCut);
    pol = sign(mean(P(indexCut)));
    u = U(indexCut);
    y = pol*detrend(abs(torqFit(indexCut)),'constant');

    if (inj2)
        for ind = 1:length(cut)
            % Sync basis function and data
            galIndex = find(floor(timSun-cut(ind)/24/3600)==0,1);
            y(ind) = y(ind) + P(ind)*injAmp*inSun(galIndex);
        end
    end

    if not(isempty(cut))   
        if (length(y)>=lenMin)
            % Sync basis function and data 
            galIndex = [];
            for cutSun = cut
                galIndex = [galIndex find(floor(timSun-cutSun/24/3600)==0,1)];
            end
            % Linear least squares fitting to basis functions and offset
%             x = [inSun(sunIndex) outSun(sunIndex) ones(length(cut),1)];
            x = [inSun(galIndex) outSun(galIndex)];
            if not(isempty(x))
                a = inv(x'*x)*x'*y';
            
                if (and(a(1)~=0,a(2)~=0))
                    cSun = [cSun a(1)];
                    sSun = [sSun a(2)];
                    uSun = [uSun std(a'*x'-y)];
                    timSunFit = [timSunFit mean(timFit(indexCut))];
                    longFit = [longFit a'*x'];
                    longTim = [longTim timFit(indexCut)];
                    longDat = [longDat y];
                    longU = [longU u];
        
                    longFit = [longFit nan];
                    longTim = [longTim nan];
                    longDat = [longDat nan];
                    longU = [longU nan];
        
                end
            end
        end
    end
end

torqSun = cSun+i*sSun;

% ampSun = sqrt(mean(cSun).^2+mean(sSun).^2);
ampSun = mean(cSun);
uncSun = std(cSun)/sqrt(length(cSun));

etaSun = ampSun/(r*m*aSun);
etaSunUnc = uncSun/(r*m*aSun);

disp(['Cosine Sun: ' num2str(mean(cSun)*1e15) ' fN m +- ' num2str(std(cSun)/sqrt(length(cSun))*1e15) ' fN m'])
disp(['Sine Sun: ' num2str(mean(sSun)*1e15) ' fN m +- ' num2str(std(sSun)/sqrt(length(sSun))*1e15) ' fN m'])
disp(['Amp Sun: ' num2str(ampSun*1e15) ' fN m +- ' num2str(uncSun*1e15) ' fN m'])
disp([' '])
disp(['Eta Sun: ' num2str(etaSun) ' +- ' num2str(etaSunUnc)])
disp([' '])

%% Galaxy Fits

cGal = [];
sGal = [];
timGalFit = [];
uGal = [];

longFit = [];
longTim = [];
longDat = [];
longU = [];

fDay = 1/24/3600;
wDay = 2*pi*fDay;
daySamples = floor(sampF/fDay);
numDaysFit = 2;
lenMin = 0.75*24*3600*TTFreq/2;

for index = 0:lenDays/numDaysFit

    indexCut = find(floor((timFit-timFit(1))/24/3600/numDaysFit)==index);
    cut = timFit(indexCut);
    pol = sign(mean(P(indexCut)));
    u = U(indexCut);
    y = pol*detrend(abs(torqFit(indexCut)),'constant');

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
            % Linear least squares fitting to basis functions and offset
%             x = [inSun(sunIndex) outSun(sunIndex) ones(length(cut),1)];
            x = [inGal(galIndex) outGal(galIndex)];
            galLongIndex = find(and(timGal>=(index*numDaysFit+timFit(1)/24/3600),...
                timGal<=((index+1)*numDaysFit)+timFit(1)/24/3600));
            xLong = [inGal(galLongIndex) outGal(galLongIndex)];
            if not(isempty(x))
                a = inv(x'*x)*x'*y';
            
                if (and(a(1)~=0,a(2)~=0))
                    cGal = [cGal a(1)];
                    sGal = [sGal a(2)];
                    uGal = [uGal std(a'*x'-y)];
                    timGalFit = [timGalFit; timGal(galLongIndex)];
                    longFit = [longFit a'*xLong'];
                    longTim = [longTim timFit(indexCut)];
                    longDat = [longDat y];
                    longU = [longU u];
        
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

torqGal = cGal+i*sGal;

ampGal = mean(cGal);
uncGal = std(cGal)/sqrt(length(cGal));

etaGalaxy = ampGal/(r*m*aGalaxy);
etaGalaxyUnc = uncGal/(r*m*aGalaxy);

disp(['Cosine Galaxy: ' num2str(mean(cGal)/(r*m)*1e15) ' fm/s^2 +- ' num2str(std(cGal)/(r*m)/sqrt(length(cGal))*1e15) ' fm/s^2'])
disp(['Sine Galaxy: ' num2str(mean(sGal)/(r*m)*1e15) ' fm/s^2 +- ' num2str(std(sGal)/(r*m)/sqrt(length(sGal))*1e15) ' fm/s^2'])
disp(['Amp Galaxy: ' num2str(ampGal/(r*m)*1e15) ' fm/s^2 +- ' num2str(uncGal/(r*m)*1e15) ' fm/s^2'])
disp([' '])
disp(['Eta Galaxy: ' num2str(etaGalaxy) ' +- ' num2str(etaGalaxyUnc)])

%% Figures
figure(4)
subplot(2,1,1)
l=errorbar(timFitin/3600/24, Cin*1e15,Uin*1e15,'.');
hold on
l=errorbar(timFit/3600/24, real(torqFit)*1e15,U*1e15,'.');
hold off
ylabel('Cosine Torque (fN m)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
grid on
subplot(2,1,2)
l=errorbar(timFitin/3600/24, Sin*1e15,Uin*1e15,'.');
hold on
ll=errorbar(timFit/3600/24, imag(torqFit)*1e15,U*1e15,'.');
hold off
ylabel('Sine Torque (fN m)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(ll,'MarkerSize',16);
% ylim([1e-18 1e-9])
% xlim([152 1.01*max(timFitin/3600/24)])
grid on

figure(5)
l=plot(real(torqFit)*1e15,imag(torqFit)*1e15,'.',real(thermCirc)*1e15,imag(thermCirc)*1e15);
ylabel('Sine Torque (fN m)','Interpreter', 'latex')
xlabel('Cosine Torque (fN m)','Interpreter', 'latex')
legend('Data','Thermal Noise','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(l,'LineWidth',1.5);
ylim([-0.2 0.2])
xlim([-0.2 0.2])
grid on

figure(6)
% l=errorbar(longTim/3600/24, longDat*1e15, longU*1e15,'.');
l=plot(longTim/3600/24, longDat*1e15, '.');
hold on
ll=plot(timGalFit,longFit*1e15);
hold off
ylabel('Torque Amplitude (fN m)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(ll,'LineWidth',2);
ylim([-0.02 0.02])
xlim([0.99*min(longTim/3600/24) 1.01*max(longTim/3600/24)])
legend('Data','Fit','Interpreter', 'latex')
grid on

% figure(7)
% subplot(2,1,1)
% l=errorbar(timGalFit/3600/24, real(torqGal)*1e15,uGal*1e15,'.');
% ylabel('In-Phase Daily (fN m)','Interpreter', 'latex')
% xlabel('Time (days)','Interpreter', 'latex')
% set(gca,'FontSize',16);
% % ylim([-3 3])
% grid on
% subplot(2,1,2)
% ll=errorbar(timGalFit/3600/24, imag(torqGal)*1e15,uGal*1e15,'.');
% ylabel('Out-of-Phase Daily (fN m)','Interpreter', 'latex')
% xlabel('Time (days)','Interpreter', 'latex')
% set(gca,'FontSize',16);
% set(l,'MarkerSize',16);
% set(ll,'MarkerSize',16);
% % ylim([-3 3])
% % xlim([1e-5 1e0])
% grid on

uncPhi = linspace(0,2*pi,100);
uncCirc = (std(cSun)/sqrt(length(cSun))*cos(thermPhi)+mean(cSun))...
        +i*(std(sSun)/sqrt(length(sSun))*sin(thermPhi)+mean(sSun));

figure(7)
% l=errorbar(real(torqSun)*1e15,imag(torqSun)*1e15,uSun*1e15,uSun*1e15,uSun*1e15,uSun*1e15,'.');
l=plot(real(torqSun)*1e15/(r*m),imag(torqSun)*1e15/(r*m),'.');
hold on
ll=plot(mean(cSun)*1e15/(r*m),mean(sSun)*1e15/(r*m),'+',...
    real(uncCirc)*1e15/(r*m),imag(uncCirc)*1e15/(r*m), 'Color', [0.8500 0.3250 0.0980]);
hold off
ylabel('Out-of-Phase Sun (fm/s$^2$)','Interpreter', 'latex')
xlabel('In-Phase Sun (fm/s$^2$)','Interpreter', 'latex')
text(-2.5e-3,11,['$\eta_\odot$ = ' num2str(etaSun,2) ' $\pm$ ' num2str(etaSunUnc,2)],'Interpreter', 'latex','FontSize',16)
% text(-2.5e-3,16e-3,['$\eta_{DM}$ = ' num2str(etaGalaxy,2) ' $\pm$ ' num2str(etaGalaxyUnc,2)],'Interpreter', 'latex','FontSize',16)
set(gca,'FontSize',16);
set(l,'MarkerSize',20);
set(ll,'LineWidth',2);
set(ll,'MarkerSize',12);
ylim([-15 15])
xlim([-15 15])
grid on

% figure(8)
% % l=errorbar(real(torqGal)*1e15,imag(torqGal)*1e15,uGal*1e15,uGal*1e15,uGal*1e15,uGal*1e15,'.');
% l=plot(real(torqGal)*1e15,imag(torqGal)*1e15,'.');
% ylabel('Out-of-Phase Galaxy (fN m)','Interpreter', 'latex')
% xlabel('In-Phase Galaxy (fN m)','Interpreter', 'latex')
% % text(-2.5e-3,18e-3,['$\eta_\odot$ = ' num2str(etaSun,2) ' $\pm$ ' num2str(etaSunUnc,2)],'Interpreter', 'latex','FontSize',16)
% text(-2.5e-3,7e-3,['$\eta_{DM}$ = ' num2str(etaGalaxy,2) ' $\pm$ ' num2str(etaGalaxyUnc,2)],'Interpreter', 'latex','FontSize',16)
% set(gca,'FontSize',16);
% set(l,'MarkerSize',20);
% ylim([-0.01 0.01])
% xlim([-0.01 0.01])
% grid on
%%
uncPhi = linspace(0,2*pi,100);
uncCirc = (std(cGal)/sqrt(length(cGal))*cos(thermPhi)+mean(cGal))...
        +i*(std(sGal)/sqrt(length(sGal))*sin(thermPhi)+mean(sGal));

figure(8)
% l=errorbar(real(torqGal)*1e15,imag(torqGal)*1e15,uGal*1e15,uGal*1e15,uGal*1e15,uGal*1e15,'.');
l=plot(real(torqGal)*1e15/(r*m),imag(torqGal)*1e15/(r*m),'.');
hold on
ll=plot(mean(cGal)*1e15/(r*m),mean(sGal)*1e15/(r*m),'+',...
    real(uncCirc)*1e15/(r*m),imag(uncCirc)*1e15/(r*m), 'Color', [0.8500 0.3250 0.0980]);
hold off
ylabel('Out-of-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
xlabel('In-Phase Acceleration (fm/s$^2$)','Interpreter', 'latex')
% text(-2.5e-3,18e-3,['$\eta_\odot$ = ' num2str(etaSun,2) ' $\pm$ ' num2str(etaSunUnc,2)],'Interpreter', 'latex','FontSize',16)
text(-2.5,12,['$\eta_{DM}$ = ' num2str(etaGalaxy,2) ' $\pm$ ' num2str(etaGalaxyUnc,2)],'Interpreter', 'latex','FontSize',16)
set(gca,'FontSize',16);
set(l,'MarkerSize',18);
set(ll,'LineWidth',2);
set(ll,'MarkerSize',12);
ylim([-15 15])
xlim([-15 15])
grid on
%%
% figure(9)
% l=plot(timGalFit/3600/24, abs(torqGal)*1e15,'.');
% ylabel('Daily Amplitude (fN m)','Interpreter', 'latex')
% xlabel('Time (days)','Interpreter', 'latex')
% set(gca,'FontSize',16);
% set(l,'MarkerSize',16);
% % ylim([1e-18 1e-9])
% % xlim([1e-5 1e0])
% grid on
% 
% mF = logspace(-5,-2);
% mA = mF*0+1e-25;

% figure(10)
% l=loglog(dmFreq,ampDM, mF, mA, '--', fDEP, aDEP, '--');
% ylabel('$g_{ULDM}$','Interpreter', 'latex')
% xlabel('Frequency (Hz)','Interpreter', 'latex')
% legend('Data', 'MICROSCOPE','DarkEP','Raw Torq Limits','Interpreter', 'latex')
% set(gca,'FontSize',16);
% set(l,'LineWidth',1.5);
% ylim([5e-27 5e-24])
% xlim([1e-5 1e-3])
% grid on
% 
% figure(11)
% l=plot(timSun, inSun, timSun, outSun);
% xlabel('Time (days)','Interpreter', 'latex')
% % legend('Data', 'MICROSCOPE','DarkEP','Raw Torq Limits','Interpreter', 'latex')
% set(gca,'FontSize',16);
% set(l,'LineWidth',1.5);
% % ylim([5e-27 5e-24])
% % xlim([1e-5 1e-3])
% grid on

%% Histrogram
if (true)
    figure(91)
    [Nin,Xin] = hist(Cin,1e2);
    [N,X] = hist(C,Xin);
    bar(Xin*1e15,(Nin))
    hold on
    bar(X*1e15,(N))
    plot(X*1e15,max(N)*exp(-(X.^2)/thermAmp^2),'LineWidth',3)
    hold off
    % xlim([-1 1])
    xlabel('Cosine Torque (fNm)','Interpreter', 'latex')
    ylabel('Number','Interpreter', 'latex')
    legend('Before Cuts','After Cuts','Thermal Noise','Interpreter', 'latex')
    set(gca,'FontSize',16);
    grid on
    
    figure(92)
    [Nin,Xin] = hist(Sin,1e2);
    [N,X] = hist(S,Xin);
    bar(Xin*1e15,(Nin))
    hold on
    bar(X*1e15,(N))
    plot(X*1e15,max(N)*exp(-(X.^2)/thermAmp^2),'LineWidth',3)
    hold off
    % xlim([-1 1])
    xlabel('Sine Torque (fNm)','Interpreter', 'latex')
    ylabel('Number','Interpreter', 'latex')
    legend('Before Cuts','After Cuts','Thermal Noise','Interpreter', 'latex')
    set(gca,'FontSize',16);
    grid on
    
    sig = thermAmp/6;
    x0 = 2.1e-17;
    k = 3;

    figure(93)
    Xg = linspace(x0,5e-17,100);
    [Nin,Xin] = hist(Uin,1e2);
    [N,X] = hist(U,Xin);
    bar(Xin,(Nin))
    hold on
    bar(X,(N))
%     plot(Xg*1e15,5e-18*max([N Nin])*((Xg-x0)/sig).^(k/2-1)./sig.*exp(-((Xg-x0)/sig)/2)/(2^(k/2)*gamma(k/2)),'LineWidth',3)
    hold off
    % xlim([-1 1])
    xlabel('Uncertainty (fNm)','Interpreter', 'latex')
    ylabel('Number','Interpreter', 'latex')
    legend('Before Cuts','After Cuts','Interpreter', 'latex')
    set(gca,'FontSize',16);
    grid on
    
end

%%

figure(94)
Xg = linspace(-1e-17,1e-17,100);
Xin = linspace(-1e-17,1e-17,15);
[Nin,Xin] = hist(real(torqGal),Xin);
[N,X] = hist(imag(torqGal),Xin);
bar(Xin*1e15,(Nin))
alpha(.5)
hold on
bar(X*1e15,(N))
alpha(.5)
plot(Xg*1e15,max([N Nin])*exp(-(Xg.^2)/(thermAmp/2)^2),'LineWidth',3)
hold off
% xlim([-1 1])
xlabel('Galaxy Torque (fNm)','Interpreter', 'latex')
ylabel('Number','Interpreter', 'latex')
legend('In-Phase','Out-of-Phase','Thermal Noise','Interpreter', 'latex')
set(gca,'FontSize',16);
grid on

%%
if(true)
    fig2=figure(6);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_FitsMultiRun.pdf','-dpdf','-r1200')
    
    fig2=figure(8);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_DailyMultiRun.pdf','-dpdf','-r1200')

    fig2=figure(7);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_SunMultiRun.pdf','-dpdf','-r1200')
end
