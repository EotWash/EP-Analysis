warning('off')
%%

injAmp = 20e-5*(r*m*aGalaxy);
inj1 = false; % Quadrature injection
inj2 = false; % Amplitude injection
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

%% Data loading

if (true)

    runs = ["run6875Fits.mat" "run6891Fits.mat" "run6893Fits.mat" "run6895Fits.mat"];
%     runs = ["run6891Fits.mat"];
    timFitin =[];
    Cin = [];
    Sin = [];
    Uin = [];
    Pin = [];
    for f=1:length(runs)
        in = load(runs(f));

        unCut = find(in.out(4,:)<prctile(in.out(4,:),80));

        timFitin = [timFitin in.out(1,unCut)];
        Cin = [Cin detrend(in.out(2,unCut))];
        Sin = [Sin detrend(in.out(3,unCut))];
        Uin = [Uin in.out(4,unCut)];
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
    rawSun=load('galVectMin.out');
    sunSampF = 1/(rawSun(2,1)-rawSun(1,1))/3600/24;
    timSun=decimate(rawSun(:,1),floor(sunSampF/sampF));
    inSun=decimate(rawSun(:,2),floor(sunSampF/sampF));
    outSun=decimate(rawSun(:,3),floor(sunSampF/sampF));
 
    if (inj1)
        for ind = 1:length(timFit)
            % Sync basis function and data
            sunIndex = find(floor(timSun-timFit(ind)/24/3600)==0,1);
            C(ind) = C(ind) + sqrt(2)/2*P(ind)*injAmp*inSun(sunIndex)+injAmp;
            S(ind) = S(ind) + sqrt(2)/2*P(ind)*injAmp*inSun(sunIndex)+injAmp;
        end
    end
end

lenDays = ceil((timFit(end)-timFit(1))/24/3600);
%%

torqFit = C+i*S;

amp = sqrt(mean(C)^2+mean(S)^2);
unc = 1/amp*sqrt(mean(C)^2*std(C)^2+mean(S)^2*std(S)^2)/sqrt(length(C));

etaEarth = amp/(r*m*aEarth);
etaEarthUnc = unc/(r*m*aEarth);

disp(['Cosine Torque: ' num2str(mean(C)*1e15) ' fN m +- ' num2str(std(C)/sqrt(length(C))*1e15) ' fN m'])
disp(['Sine Torque: ' num2str(mean(S)*1e15) ' fN m +- ' num2str(std(S)/sqrt(length(S))*1e15) ' fN m'])
disp(['Amp Torque: ' num2str(amp*1e15) ' fN m +- ' num2str(unc*1e15) ' fN m'])
disp([' '])
disp(['Eta Earth: ' num2str(etaEarth) ' +- ' num2str(etaEarthUnc)])
disp([' '])

%% Thermal Circle Calculations

thermPhi = linspace(0,2*pi,100);
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq));
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

%% Dark Matter Calcuations

amp1w = sqrt(C.^2+S.^2);

sampFDM = sampF;

dmFreq = linspace(1/24/3600,0.95*sampFDM/2,floor(sampFDM*24*3600*2));

ampDM = [];
uncDM = [];

numDaysFitDM = 1;

for indexDM = 1:length(dmFreq)
    
    fFitDM = dmFreq(indexDM);
    wFit = 2*pi*fFitDM;
    fitSamplesDM = floor(sampFDM/fFitDM);
    
    CDM = [];
    SDM = [];
    
    for index = 0:lenDays/numDaysFitDM
        
        indexCut = find(floor((timFit-timFit(1))/24/3600/numDaysFitDM)==index);

        cut = timFit(indexCut)';

        if not(isempty(cut))
            nPeriods = (cut(end)-cut(1))*fFitDM;
           
            x = [cos(wFit*cut) sin(wFit*cut)];  
    %             ones(length(cut),1)];
            
            y = detrend(amp1w(indexCut))';
        
            a = inv(x'*x)*x'*y;
            if (and(and(a(1)~=0,a(2)~=0),nPeriods>=2))
                CDM = [CDM a(1)];
                SDM = [SDM a(2)];
            end 
        end
    
    end
        
    ampDM = [ampDM; sqrt(mean(CDM)^2+mean(SDM)^2)];
    uncDM = [uncDM; 1/amp*sqrt(mean(CDM)^2*std(CDM)^2+mean(SDM)^2*std(SDM)^2)/sqrt(length(CDM))];

end

ampDM = aDM*ampDM;
uncDM = aDM*uncDM;


%% Daily Fits

cDay = [];
sDay = [];
timDay = [];
uDay = [];

longFit = [];
longTim = [];
longDat = [];
longU = [];

fDay = 1/24/3600;
wDay = 2*pi*fDay;
daySamples = floor(sampF/fDay);
lenMin = 0.75*24*3600*TTFreq/2;
numDaysFit = 1;

for index = 0:lenDays/numDaysFit

    indexCut = find(floor((timFit-timFit(1))/24/3600/numDaysFit)==index);
    cut = timFit(indexCut);
    pol = sign(mean(P(indexCut)));
    u = U(indexCut);
%     uw = 1./u.^2/sum(1./u.^2)*length(u);
    y = pol*detrend(abs(torqFit(indexCut)),'constant');

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
            % Linear least squares fitting to basis functions and offset
%             x = [inSun(sunIndex) outSun(sunIndex) ones(length(cut),1)];
            x = [inSun(sunIndex) outSun(sunIndex)];
            %     x = [cos(wDay*cut') sin(wDay*cut')...        
            %             ones(length(cut),1)];
            if not(isempty(x))
                a = inv(x'*x)*x'*y';
            
                if (and(a(1)~=0,a(2)~=0))
                    cDay = [cDay a(1)];
                    sDay = [sDay a(2)];
                    uDay = [uDay std(a'*x'-y)];
                    timDay = [timDay mean(timFit(indexCut))];
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

torqDaily = cDay+i*sDay;

ampDaily = mean(cDay);
uncDaily = std(cDay)/sqrt(length(cDay));

etaSun = ampDaily/(r*m*aSun);
etaSunUnc = uncDaily/(r*m*aSun);

etaGalaxy = ampDaily/(r*m*aGalaxy);
etaGalaxyUnc = uncDaily/(r*m*aGalaxy);

disp(['Cosine Daily: ' num2str(mean(cDay)*1e15) ' fN m +- ' num2str(std(cDay)/sqrt(length(cDay))*1e15) ' fN m'])
disp(['Sine Daily: ' num2str(mean(sDay)*1e15) ' fN m +- ' num2str(std(sDay)/sqrt(length(sDay))*1e15) ' fN m'])
disp(['Amp Daily: ' num2str(ampDaily*1e15) ' fN m +- ' num2str(uncDaily*1e15) ' fN m'])
disp([' '])
disp(['Eta Daily: ' num2str(etaSun) ' +- ' num2str(etaSunUnc)])
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
% xlim([1e-5 1e0])
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
ll=plot(longTim/3600/24,longFit*1e15);
hold off
ylabel('Torque Amplitude (fN m)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(ll,'LineWidth',2);
% ylim([-0.02 0.02])
% xlim([1e-5 1e0])
legend('Data','Fit','Interpreter', 'latex')
grid on

figure(7)
subplot(2,1,1)
l=errorbar(timDay/3600/24, real(torqDaily)*1e15,uDay*1e15,'.');
ylabel('In-Phase Daily (fN m)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
% ylim([-3 3])
grid on
subplot(2,1,2)
ll=errorbar(timDay/3600/24, imag(torqDaily)*1e15,uDay*1e15,'.');
ylabel('Out-of-Phase Daily (fN m)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(ll,'MarkerSize',16);
% ylim([-3 3])
% xlim([1e-5 1e0])
grid on

figure(8)
l=errorbar(real(torqDaily)*1e15,imag(torqDaily)*1e15,uDay*1e15,uDay*1e15,uDay*1e15,uDay*1e15,'.');
% l=plot(real(torqDaily)*1e15,imag(torqDaily)*1e15,'.');
ylabel('Out-of-Phase Daily (fN m)','Interpreter', 'latex')
xlabel('In-Phase Daily (fN m)','Interpreter', 'latex')
% text(-2.5e-3,18e-3,['$\eta_\odot$ = ' num2str(etaSun,2) ' $\pm$ ' num2str(etaSunUnc,2)],'Interpreter', 'latex','FontSize',16)
% text(-2.5e-3,16e-3,['$\eta_{DM}$ = ' num2str(etaGalaxy,2) ' $\pm$ ' num2str(etaGalaxyUnc,2)],'Interpreter', 'latex','FontSize',16)
set(gca,'FontSize',16);
set(l,'MarkerSize',20);
ylim([-0.01 0.01])
xlim([-0.01 0.01])
grid on

figure(9)
l=plot(timDay/3600/24, abs(torqDaily)*1e15,'.');
ylabel('Daily Amplitude (fN m)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
% ylim([1e-18 1e-9])
% xlim([1e-5 1e0])
grid on

mF = logspace(-5,-2);
mA = mF*0+1e-25;

figure(10)
l=loglog(dmFreq,ampDM, mF, mA, '--', fDEP, aDEP, '--');
ylabel('$g_{ULDM}$','Interpreter', 'latex')
xlabel('Frequency (Hz)','Interpreter', 'latex')
legend('Data', 'MICROSCOPE','DarkEP','Raw Torq Limits','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
ylim([5e-27 5e-24])
xlim([1e-5 1e-3])
grid on

figure(11)
l=plot(timSun, inSun, timSun, outSun);
xlabel('Time (days)','Interpreter', 'latex')
% legend('Data', 'MICROSCOPE','DarkEP','Raw Torq Limits','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
% ylim([5e-27 5e-24])
% xlim([1e-5 1e-3])
grid on

%% Histrogram

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

figure(93)
[Nin,Xin] = hist(Uin,1e2);
[N,X] = hist(U,Xin);
bar(Xin*1e15,(Nin))
hold on
bar(X*1e15,(N))
hold off
% xlim([-1 1])
xlabel('Uncertainty (fNm)','Interpreter', 'latex')
ylabel('Number','Interpreter', 'latex')
legend('Before Cuts','After Cuts','Interpreter', 'latex')
set(gca,'FontSize',16);
grid on

%%
if(false)
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

    fig2=figure(10);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'ULDM_MultiRun.pdf','-dpdf','-r1200')
end
