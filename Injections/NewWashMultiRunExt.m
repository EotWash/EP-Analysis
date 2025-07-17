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
    runs = ["run6964FitsInj.mat"];

    timFitin =[];
    Cin = [];
    Sin = [];
    Uin = [];
    Pin = [];
    Uraw = [];

    for f=1:length(runs)
        
        in = load(runs(f));        

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
        Pin = [Pin in.out(4,unCut)*0+1];
        
    end

    %% 2-sigma cuts

    % Cosine cut
    unCut = find(and(and(Cin>prctile(Cin,5),Cin<prctile(Cin,95)),...
        and(Sin>prctile(Sin,5),Sin<prctile(Sin,95))));
%     unCut = find(Cin~=0);
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
    timGal=resample(rawGal(:,1),1,ceil(galSampF/sampF));
    inGal=resample(rawGal(:,2),1,ceil(galSampF/sampF));
    outGal=detrend(resample(rawGal(:,3),1,ceil(galSampF/sampF)));
%     outGal = max(inGal)/max(outGal)*outGal;

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
%     indexCut = find(and((timFit-timFit(1))>=index*numDaysFit*24*3600,...
%         (timFit-timFit(1))<(index+1)*numDaysFit*24*3600));
%     indexCut = (index*daySamples*numDaysFit+1:(index+1)*daySamples*numDaysFit+1);

    % Cut vectors
    cut = timFit(indexCut);
    pol = sign(mean(P(indexCut)));
    u = U(indexCut);
    y = torqFit(indexCut);

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

%             % Sync basis function and data 
            galIndex = [];
            for cutGal = cut
                [m,minI] = min(abs(timGal-cutGal/24/3600));
                galIndex = [galIndex minI];
            end
%             
            % Design matrix
            x = [inGal(galIndex)+i*outGal(galIndex)];
%             x = [inGal(galIndex)+i*outGal(galIndex) (1+i)+0*galIndex'];
            
            if not(isempty(x))
               
                % Linear least squares fitting to basis functions
                a = inv(x'*x)*x'*y';

%                 figure(index+100)
%                 plot(cut,real(y))
%                 hold on
%                 plot(cut,real(a)'*real(x'))
%                 hold off
            
                % Append valid data points
                if ((a(1)~=0))
                    cGal = [cGal real(a(1))];
                    sGal = [sGal imag(a(1))];
                    uGal = [uGal std(a'*x'-y)];
                    timGalFit = [timGalFit; timGal(galIndex)];
                    longFit = [longFit (real(a'*x'))];
                    longTim = [longTim timFit(indexCut)];
%                     longDat = [longDat (real(y)+imag(y))/2];
                    longDat = [longDat real(y)];
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


%% Save plots


