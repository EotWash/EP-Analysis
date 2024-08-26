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

%%

injAmp = 1e-25/aDM;
wInj = 2*pi*6e-5;
inj1 = false; % Quadrature injection
inj2 = false; % Amplitude injection

%% Data loading

if (true)

    runs = ["run6875Fits.mat" "run6891Fits.mat" "run6893Fits.mat" ...
        "run6895Fits.mat" "run6896Fits.mat" "run6897Fits.mat" "run6900Fits.mat"];
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
    rawX=load('xVectMin.out');
    xSampF = 1/(rawX(2,1)-rawX(1,1))/3600/24;
    timX=decimate(rawX(:,1),floor(xSampF/sampF));
    inX=decimate(rawX(:,2),floor(xSampF/sampF));
    outX=decimate(rawX(:,3),floor(xSampF/sampF));

    % Loading in basis funtions outputted from sunVect.py
    rawY=load('yVectMin.out');
    ySampF = 1/(rawY(2,1)-rawY(1,1))/3600/24;
    timY=decimate(rawY(:,1),floor(ySampF/sampF));
    inY=decimate(rawY(:,2),floor(ySampF/sampF));
    outY=decimate(rawY(:,3),floor(ySampF/sampF));

    % Loading in basis funtions outputted from sunVect.py
    rawZ=load('zVectMin.out');
    zSampF = 1/(rawZ(2,1)-rawZ(1,1))/3600/24;
    timZ=decimate(rawZ(:,1),floor(zSampF/sampF));
    inZ=decimate(rawZ(:,2),floor(zSampF/sampF));
    outZ=decimate(rawZ(:,3),floor(zSampF/sampF));
 
    if (inj1)
        for ind = 1:length(timFit)
            % Sync basis function and data
            sunIndex = find(floor(timX-timFit(ind)/24/3600)==0,1);
            C(ind) = C(ind) + sqrt(2)/2*P(ind)*injAmp*inX(sunIndex)*cos(wInj*timX(sunIndex))+injAmp;
            S(ind) = S(ind) + sqrt(2)/2*P(ind)*injAmp*inX(sunIndex)*cos(wInj*timX(sunIndex))+injAmp;
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

%% Thermal Circle Calculations

thermPhi = linspace(0,2*pi,100);
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq));
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

%% Dark Matter Calcuations

amp1w = sqrt(C.^2+S.^2);

sampFDM = sampF;

dmFreq = linspace(1/24/3600,0.95*sampFDM/2,floor(sampFDM*24*3600*2));

ampDMX = [];
uncDMX = [];

ampDMY = [];
uncDMY = [];

ampDMZ = [];
uncDMZ = [];

numDaysFitDM = 1;

for indexDM = 1:length(dmFreq)
    
    fFitDM = dmFreq(indexDM);
    wFit = 2*pi*fFitDM;
    fitSamplesDM = floor(sampFDM/fFitDM);
    
    CDMX = [];
    SDMX = [];

    CDMY = [];
    SDMY = [];

    CDMZ = [];
    SDMZ = [];
    
    for index = 0:lenDays/numDaysFitDM
        
        indexCut = find(floor((timFit-timFit(1))/24/3600/numDaysFitDM)==index);

        cut = timFit(indexCut)';

        y = detrend(amp1w(indexCut))';        

        if not(isempty(cut))
            nPeriods = (cut(end)-cut(1))*fFitDM;

            basisIndex = [];
            for cutX = cut'
                basisIndex = [basisIndex find(floor(timX-cutX/24/3600)==0,1)];
            end

            if (inj2)
                y = y + injAmp*inX(basisIndex).*cos(wInj*cut);
            end
           
%             x = [inX(basisIndex).*cos(wFit*cut) inX(basisIndex).*sin(wFit*cut)...
%                 inY(basisIndex).*cos(wFit*cut) inY(basisIndex).*sin(wFit*cut)...
%                 inZ(basisIndex).*cos(wFit*cut) inZ(basisIndex).*sin(wFit*cut)];  

            x = [inX(basisIndex).*cos(wFit*cut) inX(basisIndex).*sin(wFit*cut)...
                inY(basisIndex).*cos(wFit*cut) inY(basisIndex).*sin(wFit*cut)]; 

                    
            a = inv(x'*x)*x'*y;
            if (and(and(a(1)~=0,a(2)~=0),nPeriods>=2))
                CDMX = [CDMX a(1)];
                SDMX = [SDMX a(2)];

                CDMY = [CDMY a(3)];
                SDMY = [SDMY a(4)];

%                 CDMZ = [CDMZ a(5)];
%                 SDMZ = [SDMZ a(6)];
            end 
        end
    
    end
        
    ampDMX = [ampDMX; sqrt(mean(CDMX)^2+mean(SDMX)^2)];
    uncDMX = [uncDMX; 1/amp*sqrt(mean(CDMX)^2*std(CDMX)^2+mean(SDMX)^2*std(SDMX)^2)/sqrt(length(CDMX))];

    ampDMY = [ampDMY; sqrt(mean(CDMY)^2+mean(SDMY)^2)];
    uncDMY = [uncDMY; 1/amp*sqrt(mean(CDMY)^2*std(CDMY)^2+mean(SDMY)^2*std(SDMY)^2)/sqrt(length(CDMY))];
    

    ampDMZ = [ampDMZ; sqrt(mean(CDMZ)^2+mean(SDMZ)^2)];
    uncDMZ = [uncDMZ; 1/amp*sqrt(mean(CDMZ)^2*std(CDMZ)^2+mean(SDMZ)^2*std(SDMZ)^2)/sqrt(length(CDMZ))];
end

ampDMX = aDM*ampDMX;
uncDMX = aDM*uncDMX;

ampDMY = aDM*ampDMY;
uncDMY = aDM*uncDMY;

ampDMZ = aDM*ampDMZ;
uncDMZ = aDM*uncDMZ;

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

mF = logspace(-5,-2);
mA = mF*0+1e-25;

figure(10)
l=loglog(dmFreq,ampDMX, dmFreq,ampDMY, dmFreq,ampDMZ, mF, mA, '--', fDEP, aDEP, '--');
ylabel('$g_{ULDM}$','Interpreter', 'latex')
xlabel('Frequency (Hz)','Interpreter', 'latex')
legend('Data X', 'Data Y', 'Data Z', 'MICROSCOPE','DarkEP','Raw Torq Limits','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
ylim([5e-27 5e-24])
xlim([1e-5 1e-3])
grid on



%%
if(false)

    fig2=figure(10);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'ULDM_MultiRun.pdf','-dpdf','-r1200')
end
