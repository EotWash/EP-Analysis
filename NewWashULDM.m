warning('off')

%%

w0 = 2*pi*6.8567e-4; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
Q = 1.1e5; % Quality factor
kappa = I*w0^2; % Spring constance (N m/rad)
kb = 1.38064852e-23; % Boltzmann's constant (J/K)
T = 293; % Temperature (K)
thetaCalib = 3/300/8; % Autocollimator calibration (rad/(Diff/Sum))
m = 38.72e-3/2; % Mass (kg)
r = 3.77e-2/2; % Lever-arm (m)

f2M = 4.1e-19/1e-4;

TTFreq = 0.457120e-3; % Turn table frequency (Hz)
wTT = 2*pi*TTFreq; % Turn table frequency (rad*Hz)

aDM = 1e-25/9.73e-19; % Torque to g_dm conversion for Be-Al

% Thermal noise
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq)); 

%% Injection Controls

injAmp = 1e-23/aDM;
wInj = 2*pi*2e-19/f2M;
inj1 = false; % Quadrature injection
inj2 = true; % Amplitude injection

%% Data loading

if (true)

    % Runs to load. Once per turntable cosine amplitude, sine amplitude, 
    % and misfit are calculated in NewWashAnalysis.m then loaded here
    runs = ["run6891Fits.mat" "run6893Fits.mat" ...
        "run6895Fits.mat" "run6896Fits.mat" "run6897Fits.mat" "run6900Fits.mat" ...
        "run6903Fits.mat" "run6904Fits.mat" "run6905Fits.mat" "run6923Fits.mat" ...
        "run6925Fits.mat" "run6926Fits.mat" "run6927Fits.mat" "run6930Fits.mat"...
        "run6931Fits.mat" "run6936Fits.mat" "run6939Fits.mat" "run6949Fits.mat"...
        "run6950Fits.mat"];

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
        if f<10
            timFitin = [timFitin mod(in.out(1,unCut),31556926)];
        else
            timFitin = [timFitin mod(in.out(1,unCut),31556926)+31556926];
        end
        Cin = [Cin detrend(in.out(2,unCut))];
        Sin = [Sin detrend(in.out(3,unCut))];
        Uin = [Uin in.out(4,unCut).^2/(thermAmp^2)/sqrt(length(in.out(4,:)))];
        Uraw = [Uraw in.out(4,:).^2/(thermAmp^2)/sqrt(length(in.out(4,:)))];

        % Pendulum flips
        if or(f<3,and(f>13,f<20))
            Pin = [Pin in.out(4,unCut)*0+1];
        else 
            Pin = [Pin in.out(4,unCut)*0-1];
        end
        if f==2
            timFlip = mod(in.out(1,end),31556926);
        end
    end
    
    % Moved time zero to midnight Jan. 1, 2024
%     timFitin = mod(timFitin,31556926);

    % Load previouse limits from DarkEP
    inDEP = load('DarkEPData.csv');
    fDEP = inDEP(:,1);
    aDEP = inDEP(:,2);

    % Load previouse limits from LISA Pathfinder
    inLISA = load('LISA-PF.csv');
    fLISA = inLISA(:,1);
    aLISA = inLISA(:,2)/0.428;     

    %% 2-sigma cuts

    % Cosine cut
    unCut = find(and(and(Cin>prctile(Cin,5),Cin<prctile(Cin,95)),...
        and(Sin>prctile(Sin,5),Sin<prctile(Sin,95))));
    C = Cin(unCut);
    S = Sin(unCut);
    U = Uin(unCut);
    P = Pin(unCut);
    timFit = timFitin(unCut);

    % Sine cut
%     unCut = find(and(S>prctile(S,5),S<prctile(S,95)));
%     C = C(unCut);
%     S = S(unCut);
%     U = U(unCut);
%     P = P(unCut);
%     timFit = timFit(unCut);
    
    % Sampling frequency
    sampF = 1/(timFit(5)-timFit(4));

    % Loading in X basis funtions outputted from dmVect.py
    rawX=load('xVectMin.out');
    xSampF = 1/(rawX(2,1)-rawX(1,1))/3600/24;
    timX=decimate(rawX(:,1),floor(xSampF/sampF));
    inX=decimate(rawX(:,2),floor(xSampF/sampF));
    outX=decimate(rawX(:,3),floor(xSampF/sampF));

    timX = timX - (timX>307.042)/24 - (timX<69.082)/24 - (timX<433.041)/24 - (timX>671.041)/24;

    timX = [timX; timX+365.2422];
    inX = [inX; inX];
    outX = [outX; outX];

    % Loading in Y basis funtions outputted from dmVect.py
    rawY=load('yVectMin.out');
    ySampF = 1/(rawY(2,1)-rawY(1,1))/3600/24;
    timY=decimate(rawY(:,1),floor(ySampF/sampF));
    inY=decimate(rawY(:,2),floor(ySampF/sampF));
    outY=decimate(rawY(:,3),floor(ySampF/sampF));

    timY = timY - (timY>307.042)/24 - (timY<69.082)/24 - (timY<433.041)/24 - (timY>671.041)/24;

    timY = [timY; timY+365.2422];
    inY = [inY; inY];
    outY = [outY; outY];

    % Loading in Z basis funtions outputted from dmVect.py
    rawZ=load('zVectMin.out');
    zSampF = 1/(rawZ(2,1)-rawZ(1,1))/3600/24;
    timZ=decimate(rawZ(:,1),floor(zSampF/sampF));
    inZ=decimate(rawZ(:,2),floor(zSampF/sampF));
    outZ=decimate(rawZ(:,3),floor(zSampF/sampF));

    timZ = timZ - (timZ>307.042)/24 - (timZ<69.082)/24 - (timZ<433.041)/24 - (timZ>671.041)/24;

    timZ = [timZ; timZ+365.2422];
    inZ = [inZ; inZ];
    outZ = [outZ; outZ];
 
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

% Minimum number of periods in each cut
periodMin = 0;

% Number of periods to fit
fitPeriods = 7;

% Thermal noise circle
thermPhi = linspace(0,2*pi,100); 
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

%% Fits

% Amplitude of 1 omega torque
% amp1w = sqrt(C.^2+S.^2);
% amp1w = thermAmp/100*randn(size(C));
% timFit = linspace(min(timFit),max(timFit),length(C));

% Dark matter search frequencies
dmFreq = linspace(1/24/3600/lenDays*2,sampF*2,300)';

% Create vectors
ampDMX = [];
uncDMX = [];
ampDMY = [];
uncDMY = [];
ampDMZ = [];
uncDMZ = [];
longC = [];
longS = [];

% Fit for each frequency in vector
for indexDM = 1:length(dmFreq)
    
    longAmp = [];
    longTim = [];

    % Fit parameters
    fFitDM = dmFreq(indexDM);
    wFit = 2*pi*fFitDM;
    fitSamplesDM = floor(sampF/dmFreq(1)/fitPeriods);
    
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
            pol = sign(mean(P(indexCut)));
            y = pol*detrend(abs(torqFit(indexCut)),'constant');        
    
            if not(isempty(cut))                
                
                % Sync basis function and data
                basisIndex = [];
                for cutX = cut'
                    basisIndex = [basisIndex find(floor(timX-cutX/24/3600)==0,1)];
                end
                
                % Amplitude injection
                if (inj2)
                    y = (y' + injAmp*inX(basisIndex).*cos(wInj*cut+pi/4))';
                end
                
                % Design matrix
                x = [inX(basisIndex).*cos(wFit*cut) inX(basisIndex).*sin(wFit*cut)...
                    inY(basisIndex).*cos(wFit*cut) inY(basisIndex).*sin(wFit*cut)...
                    inZ(basisIndex).*cos(wFit*cut) inZ(basisIndex).*sin(wFit*cut)]; 
    
                % Linear least squares fitting to basis functions  
                a = inv(x'*x)*x'*y';
                
                % Calculate number of periods in cut
                nPeriods = (cut(end)-cut(1))*fFitDM;
                
                % Append valid data points
                if (and(and(a(1)~=0,a(2)~=0),nPeriods>=periodMin))
                    CDMX = [CDMX a(1)];
                    SDMX = [SDMX a(2)];     
                    CDMY = [CDMY a(3)];
                    SDMY = [SDMY a(4)];     
                    CDMZ = [CDMZ a(5)];
                    SDMZ = [SDMZ a(6)];
                end 
            end
        end
        longAmp = [longAmp y];
        longTim = [longTim cut'];
    end
    
    % Calculate amplitude and uncertainty for each direction
    ampDMX = [ampDMX; sqrt(mean(CDMX)^2+mean(SDMX)^2)];
    uncDMX = [uncDMX; 1/sqrt(mean(CDMX)^2+mean(SDMX)^2)*sqrt(mean(CDMX)^2*std(CDMX)^2+mean(SDMX)^2*std(SDMX)^2)/sqrt(length(CDMX))];
    ampDMY = [ampDMY; sqrt(mean(CDMY)^2+mean(SDMY)^2)];
    uncDMY = [uncDMY; 1/sqrt(mean(CDMY)^2+mean(SDMY)^2)*sqrt(mean(CDMY)^2*std(CDMY)^2+mean(SDMY)^2*std(SDMY)^2)/sqrt(length(CDMY))];
    ampDMZ = [ampDMZ; sqrt(mean(CDMZ)^2+mean(SDMZ)^2)];
    uncDMZ = [uncDMZ; 1/sqrt(mean(CDMZ)^2+mean(SDMZ)^2)*sqrt(mean(CDMZ)^2*std(CDMZ)^2+mean(SDMZ)^2*std(SDMZ)^2)/sqrt(length(CDMZ))];

    longC = [longC CDMX CDMY CDMZ];
    longS = [longS SDMX SDMY SDMZ];   
end

% Convert to g_B-L units
ampDMX = aDM*ampDMX;
uncDMX = aDM*uncDMX*2;
ampDMY = aDM*ampDMY;
uncDMY = aDM*uncDMY*2;
ampDMZ = aDM*ampDMZ;
uncDMZ = aDM*uncDMZ*2;


%% Figures

% Time Series
Rat = 9;
figure(1)
subplot(1,Rat,[1 Rat-1])
% ll=plot(timFit/3600/24, detrend(abs(torqFit),'constant')*1e18,'.',...
ll=plot(longTim/3600/24, longAmp*1e18,...
    [213 213],[-40 40],'k--', [420 420],[-40 40],'k--');
% hold on
% patch([275 384 384 275], [-20 -20 20 20], [.5 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5)
% text(282, 2, 'Hardware Failures','Interpreter', 'latex','FontSize',14)
% text(195, 18, '0$^\circ$','Interpreter', 'latex','FontSize',16)
% text(235, 18, '180$^\circ$','Interpreter', 'latex','FontSize',16)
% text(390, 18, '180$^\circ$','Interpreter', 'latex','FontSize',16)
% text(455, 18, '0$^\circ$','Interpreter', 'latex','FontSize',16)
% hold off
ylabel('Torque (aN m)','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
set(gca,'FontSize',18);
set(ll,'MarkerSize',16);
set(ll,'LineWidth',2);
ylim([-20 20])
xlim([0.99*min(timFit/3600/24) 1.01*max(timFit/3600/24)])
grid on

subplot(1,Rat,Rat)
[n,x] = hist(longAmp*1e18);
barh(x,n,1);
% ylim([-30 30])
ylim([-20 20])
% xlim([-7 180])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
set(gca,'XGrid','off','YGrid','on')
%%
% Limits mass x-axis
figure(2)
l=loglog(dmFreq,ampDMX, dmFreq, ampDMY, dmFreq,ampDMZ);
ylabel('$g_{B-L}/\sqrt{\hbar c}$','Interpreter', 'latex')
xlabel('Frequency (Hz)','Interpreter', 'latex')
legend('X-Limits', 'Y-Limits', 'Z-Limits', 'MICROSCOPE','Shaw et. al.','LISA Pathfinder','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
% ylim([9e-27 1.1e-24])
grid on
%%
% Limits mass x-axis
mF = logspace(-7.5,-2.5);
mA = mF*0+2.753e-25;
proj = min([mF*0+7e-27; sqrt((4e-27*sqrt(1e-3)*1./sqrt(mF)).^2+ (1e-26/1e-1*mF).^2)]);
dmAmp = sqrt(uncDMX.^2 + uncDMY.^2 + uncDMZ.^2)/sqrt(3);
dmAmpA = sqrt(ampDMX.^2 + ampDMY.^2 + ampDMZ.^2)/sqrt(3);
dmFreqP = dmFreq;

dmIndex = find(and(not(isnan(dmAmp)), f2M*dmFreqP>2.2e-21));
dmAmpPlot = dmAmp(dmIndex)';
dmFreqPlot = dmFreqP(dmIndex)';

lIndex = find(f2M*fLISA<2.2e-17);
aLISAPlot = aLISA(lIndex)';
fLISAPlot = fLISA(lIndex)';

%%
figure(3)
l=loglog(f2M*dmFreqP, dmAmp, f2M*dmFreqP, movmean(dmAmp,500),...    
    f2M*mF, mA, '--', f2M*fDEP, aDEP, '--',[2.2e-21 2.2e-21],[1e-27 1e-24] ,'--', f2M*fLISA,aLISA,'--', ...
    f2M*[1/24/3600 1/24/3600], [1e-28 1e-25],'k-.');
% hold on
% patch([f2M*mF fliplr(f2M*mF)], [mA 1e-24*ones(size(mA))], [.7 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5) 
% patch([f2M*fLISAPlot fliplr(f2M*fLISAPlot)], [aLISAPlot mean(mA)*ones(size(aLISAPlot))], [.7 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5)
% patch([f2M*dmFreqPlot fliplr(f2M*dmFreqPlot)], [dmAmpPlot mean(mA)*ones(size(dmAmpPlot))], [.7 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5)
% patch([[3.58e-22 1.79e-21 2.2e-21] [2.2e-21 1.79e-21 3.58e-22]], [2.753e-25 2.753e-25 2.753e-25 1e-27 1e-27 1e-27], [.7 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5)
% patch([[2.2e-21 2.2442e-21] [2.2442e-21 2.2e-21]], [2.753e-25 2.753e-25 9.23e-26 9.23e-26], [.7 .7 .7], 'LineStyle', 'none', 'FaceAlpha', 0.5)
text(3.8e-20, 6e-27, 'Daily Frequency','Interpreter', 'latex','FontSize',16,'Rotation',90)
% text(3.8e-20/2, 4e-27, 'Half Daily Frequency','Interpreter', 'latex','FontSize',16,'Rotation',90)
% hold off
ylabel('$g_{B-L}/\sqrt{\hbar c}$','Interpreter', 'latex')
xlabel('Mass (eV)','Interpreter', 'latex')
legend('Amplitude Limits', 'Smoothed Limits', 'MICROSCOPE','Shaw et. al.','Zimmermann et. al.','LISA Pathfinder','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
% ylim([5e-27 1e-24])
xlim([4e-22 1e-17])
grid on

figure(4)
[NC,XC] = hist(aDM*longC,50);
[NS,XS] = hist(aDM*longS,XC);
bar(XC,(NC),'FaceAlpha',0.5)
hold on
bar(XS,(NS),'FaceAlpha',0.5)
hold off
xlim([-4 4]*1e-25)
xlabel('$g_{B-L}/\sqrt{\hbar c}$','Interpreter', 'latex')
ylabel('Number','Interpreter', 'latex')
legend('Cosine','Sine','Thermal Noise','Interpreter', 'latex')
set(gca,'FontSize',16);
grid on
%%
% Remove insanely large chi-squareds
Uraw = Uraw(find(Uraw<10));

figure(5)
[NUr,XU] = hist(Uraw,50);
[NU,XU] = hist(U,XU);
bar(XU,(NUr),'FaceAlpha',0.5)
hold on
bar(XU,(NU),'FaceAlpha',0.5)
plot([4 4],[0 500],'k--','LineWidth',1.5)
hold off
% xlim([-4 4]*1e-25)
xlabel('$\chi^2$ Relative to Thermal Noise','Interpreter', 'latex')
ylabel('Number','Interpreter', 'latex')
legend('Before Cuts','After Cuts','Threshold','Interpreter', 'latex')
set(gca,'FontSize',16);
grid on

%% Save plots

if(false)

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

    fig2=figure(5);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'ULDM_Hist.pdf','-dpdf','-r1200')

end
