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

TTFreq = 0.457120e-3; % Turn table frequency (Hz)
wTT = 2*pi*TTFreq; % Turn table frequency (rad*Hz)

aDM = 1e-25/9.73e-19; % Torque to g_dm conversion for Be-Al

% Thermal noise
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq)); 

%% Injection Controls

injAmp = 1e-24/aDM;
wInj = 2*pi*6e-5;
inj1 = false; % Quadrature injection
inj2 = false; % Amplitude injection

%% Data loading

if (true)

    % Runs to load. Once per turntable cosine amplitude, sine amplitude, 
    % and misfit are calculated in NewWashAnalysis.m then loaded here
    runs = ["run6875Fits.mat" "run6891Fits.mat" "run6893Fits.mat" ...
        "run6895Fits.mat" "run6896Fits.mat" "run6897Fits.mat" "run6900Fits.mat" ...
        "run6903Fits.mat" "run6904Fits.mat" "run6905Fits.mat"];

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
        Uin = [Uin in.out(4,unCut)];

        % Pendulum flip after run6895
        if f<4
            Pin = [Pin in.out(4,unCut)*0+1];
        else 
            Pin = [Pin in.out(4,unCut)*0-1];
        end
    end
    
    % Moved time zero to midnight Jan. 1, 2024
    timFitin = mod(timFitin,31556926);

    % Load previouse limits from DarkEP
    inDEP = load('DarkEPData.csv');
    fDEP = inDEP(:,1);
    aDEP = inDEP(:,2);

    % Load previouse limits from LISA Pathfinder
    inLISA = load('LISA-PF.csv');
    fLISA = inLISA(:,1);
    aLISA = inLISA(:,2);     

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

    % Loading in X basis funtions outputted from dmVect.py
    rawX=load('xVectMin.out');
    xSampF = 1/(rawX(2,1)-rawX(1,1))/3600/24;
    timX=decimate(rawX(:,1),floor(xSampF/sampF));
    inX=decimate(rawX(:,2),floor(xSampF/sampF));
    outX=decimate(rawX(:,3),floor(xSampF/sampF));

    % Loading in Y basis funtions outputted from dmVect.py
    rawY=load('yVectMin.out');
    ySampF = 1/(rawY(2,1)-rawY(1,1))/3600/24;
    timY=decimate(rawY(:,1),floor(ySampF/sampF));
    inY=decimate(rawY(:,2),floor(ySampF/sampF));
    outY=decimate(rawY(:,3),floor(ySampF/sampF));

    % Loading in Z basis funtions outputted from dmVect.py
    rawZ=load('zVectMin.out');
    zSampF = 1/(rawZ(2,1)-rawZ(1,1))/3600/24;
    timZ=decimate(rawZ(:,1),floor(zSampF/sampF));
    inZ=decimate(rawZ(:,2),floor(zSampF/sampF));
    outZ=decimate(rawZ(:,3),floor(zSampF/sampF));
 
    % Quadrature injection
    if (inj1)
        for ind = 1:length(timFit)
            % Sync basis function and data
            sunIndex = find(floor(timX-timFit(ind)/24/3600)==0,1);
            C(ind) = C(ind) + sqrt(2)/2*P(ind)*injAmp*inX(sunIndex)*cos(wInj*timX(sunIndex))+injAmp;
            S(ind) = S(ind) + sqrt(2)/2*P(ind)*injAmp*inX(sunIndex)*cos(wInj*timX(sunIndex))+injAmp;
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

% Thermal noise circle
thermPhi = linspace(0,2*pi,100); 
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

%% Fits

% Amplitude of 1 omega torque
amp1w = sqrt(C.^2+S.^2);

% Dark matter search frequencies
dmFreq = linspace(1/24/3600/lenDays*2,0.85*sampF/2,floor(0.85*sampF*24*3600*lenDays/2))';

% Create vectors
ampDMX = [];
uncDMX = [];
ampDMY = [];
uncDMY = [];
ampDMZ = [];
uncDMZ = [];

% Fit for each frequency in vector
for indexDM = 1:length(dmFreq)
    
    % Fit parameters
    fFitDM = dmFreq(indexDM);
    wFit = 2*pi*fFitDM;
    fitSamplesDM = floor(sampF/fFitDM*(fFitDM*16/(dmFreq(end)-dmFreq(1))+1));
    
    % Create vectors
    CDMX = [];
    SDMX = [];
    CDMY = [];
    SDMY = [];
    CDMZ = [];
    SDMZ = [];
    
    for index = 0:floor(length(timFit)/fitSamplesDM)
        
        % Cut indices
        indexCut = (index*fitSamplesDM+1:(index+1)*fitSamplesDM+1);
        
        % End of vector check
        if indexCut(end)<length(timFit)
            
            % Cut vectors
            cut = timFit(indexCut)';    
            y = P(indexCut).*detrend(amp1w(indexCut))';        
    
            if not(isempty(cut))                
                
                % Sync basis function and data
                basisIndex = [];
                for cutX = cut'
                    basisIndex = [basisIndex find(floor(timX-cutX/24/3600)==0,1)];
                end
                
                % Amplitude injection
                if (inj2)
                    y = y + injAmp*inX(basisIndex).*cos(wInj*cut+pi/4);
                end
                
                % Design matrix
                x = [inX(basisIndex).*cos(wFit*cut) inX(basisIndex).*sin(wFit*cut)...
                    inY(basisIndex).*cos(wFit*cut) inY(basisIndex).*sin(wFit*cut)...
                    0.7385.*cos(wFit*cut) 0.7385.*sin(wFit*cut)]; 
    
                % Linear least squares fitting to basis functions  
                a = inv(x'*x)*x'*y;
                
                % Calculate number of periods in cut
                nPeriods = (cut(end)-cut(1))*fFitDM;
                
                % Append valid data points
                if (and(and(a(1)~=0,a(2)~=0),nPeriods>=0))
                    CDMX = [CDMX a(1)];
                    SDMX = [SDMX a(2)];     
                    CDMY = [CDMY a(3)];
                    SDMY = [SDMY a(4)];     
                    CDMZ = [CDMZ a(5)];
                    SDMZ = [SDMZ a(6)];
                end 
            end
        end
    
    end
    
    % Calculate amplitude and uncertainty for each direction
    ampDMX = [ampDMX; sqrt(mean(CDMX)^2+mean(SDMX)^2)];
    uncDMX = [uncDMX; 1/sqrt(mean(CDMX)^2+mean(SDMX)^2)*sqrt(mean(CDMX)^2*std(CDMX)^2+mean(SDMX)^2*std(SDMX)^2)/sqrt(length(CDMX))];
    ampDMY = [ampDMY; sqrt(mean(CDMY)^2+mean(SDMY)^2)];
    uncDMY = [uncDMY; 1/sqrt(mean(CDMY)^2+mean(SDMY)^2)*sqrt(mean(CDMY)^2*std(CDMY)^2+mean(SDMY)^2*std(SDMY)^2)/sqrt(length(CDMY))];
    ampDMZ = [ampDMZ; sqrt(mean(CDMZ)^2+mean(SDMZ)^2)];
    uncDMZ = [uncDMZ; 1/sqrt(mean(CDMZ)^2+mean(SDMZ)^2)*sqrt(mean(CDMZ)^2*std(CDMZ)^2+mean(SDMZ)^2*std(SDMZ)^2)/sqrt(length(CDMZ))];

end

% Convert to g_B-L units
ampDMX = aDM*ampDMX;
uncDMX = aDM*uncDMX;
ampDMY = aDM*ampDMY;
uncDMY = aDM*uncDMY;
ampDMZ = aDM*ampDMZ;
uncDMZ = aDM*uncDMZ;


%% Figures

% Time Series
figure(1)
% l=plot(timFitin/3600/24, sqrt(Cin.^2+Sin.^2)*1e15,'.');
% hold on
ll=plot(timFit/3600/24, detrend(abs(torqFit))*1e15,'.');
% hold off
ylabel('Torque (fN m)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
% legend('Before Cuts','After Cuts','Interpreter', 'latex')
set(gca,'FontSize',16);
set(ll,'MarkerSize',16);
grid on

% Limits mass x-axis
figure(2)
l=loglog(dmFreq,uncDMX, dmFreq, uncDMY, dmFreq,uncDMZ);
ylabel('$g_{B-L}$','Interpreter', 'latex')
xlabel('Frequency (Hz)','Interpreter', 'latex')
legend('X-Limits', 'Y-Limits', 'Z-Limits', 'MICROSCOPE','Shaw et. al.','LISA Pathfinder','Projected Improvements','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
ylim([4e-28 1e-24])
grid on

% Limits mass x-axis
mF = logspace(-7,-1.5);
mA = mF*0+1e-25;
proj = min([mF*0+7e-27; sqrt((4e-27*sqrt(1e-3)*1./sqrt(mF)).^2+ (1e-26/1e-1*mF).^2)]);
dmAmp = sqrt(uncDMX.^2 + uncDMY.^2 + uncDMZ.^2)/sqrt(3);
f2M = 4.1e-19/1e-4;

dmIndex = find(and(not(isnan(dmAmp)), dmAmp>0));
dmAmpPlot = dmAmp(dmIndex)';
dmFreqPlot = dmFreq(dmIndex)';

lIndex = find(f2M*fLISA<2.2e-17);
aLISAPlot = aLISA(lIndex)';
fLISAPlot = fLISA(lIndex)';

figure(3)
l=loglog(f2M*dmFreq, dmAmp,...    
    f2M*mF, mA, '--', f2M*fDEP, aDEP, '--', f2M*fLISA,aLISA,'--');%,f2M*mF, proj,'-.');%, f2M*[TTFreq/2 TTFreq/2], [1e-28 1e-23],'k-.');
hold on
patch([f2M*mF fliplr(f2M*mF)], [mA 1e-24*ones(size(mA))], [.7 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5) 
patch([f2M*fLISAPlot fliplr(f2M*fLISAPlot)], [aLISAPlot 1e-25*ones(size(aLISAPlot))], [.7 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5)
patch([f2M*dmFreqPlot fliplr(f2M*dmFreqPlot)], [dmAmpPlot 1e-25*ones(size(dmAmpPlot))], [.7 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5)
hold off
ylabel('$g_{B-L}$','Interpreter', 'latex')
xlabel('Mass (eV)','Interpreter', 'latex')
legend('Amplitude Limits', 'MICROSCOPE','Shaw et. al.','LISA Pathfinder','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
ylim([4e-28 1e-24])
xlim(f2M*[min(mF) max(mF)])
grid on

%% Save plots

if(true)

    fig2=figure(2);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'ULDM_Limits.pdf','-dpdf','-r1200')

    fig2=figure(3);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'ULDM_AmpLimits.pdf','-dpdf','-r1200')

    fig2=figure(1);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'ULDM_TimeSeries.pdf','-dpdf','-r1200')

end
