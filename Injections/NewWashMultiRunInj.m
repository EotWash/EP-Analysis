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
    outGal=(decimate(rawGal(:,3),floor(galSampF/sampF)));

    timGal = timGal - (timGal>307.042)/24 - (timGal<69.082)/24 - (timGal<433.041)/24 - (timGal>671.041)/24;

end

% Length of days
lenDays = ceil((timFit(end)-timFit(1))/24/3600);

% Calculate complex torque amplitude
torqFit = P.*(C+i*S);

% %% Injection Controls
injEx = [];
injUnc = [];

inj = linspace(-10,10,21)*1e-5;

for injIndex = 1:length(inj)

    injAmp = inj(injIndex)*(r*m*aGalaxy);
    
    % Quadrature injection    
%     torqFit = injAmp*inGal';
%     timFit = timGal'*24*3600;
%     P = ones(length(inGal),1);
%     U = 0*P';
    
    
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
    
        % Sync basis function and data 
        galIndex = [];
        for cutGal = cut
            [m,minI] = min(abs(timGal-cutGal/24/3600));
            galIndex = [galIndex minI];
        end

        % Amplitude injection
        if (true)
            y = y + injAmp.*(inGal(galIndex)+i*outGal(galIndex))';
        end
    
        if not(isempty(cut))   
            if (length(y)>=lenMin)    
                
                % Design matrix
                x = [inGal(galIndex)+i*outGal(galIndex)];
                
                if not(isempty(x))
    
                    % Linear least squares fitting to basis functions
                    a = inv(x'*x)*x'*y';
                
                    % Append valid data points
                    if ((a(1)~=0))
                        cGal = [cGal real(a(1))];
                        sGal = [sGal imag(a(1))];
                        uGal = [uGal std(a'*x'-y)];                       
            
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
        
    %%
    
    % Eotvos parameters
    etaGalaxy = ampGal/(r*m*aGalaxy);
    etaGalaxyUnc = uncGal/(r*m*aGalaxy);

    injEx = [injEx etaGalaxy];
    injUnc = [injUnc etaGalaxyUnc];

end

%% Figures

figure(1)
ll = errorbar(inj*1e5, injEx*1e5, injUnc*1e5,'.');
hold on
l = plot(inj*1e5,inj*1e5+5.6e-1);
hold off
xlabel('$\eta \times 10^{-5}$ Injected','Interpreter', 'latex')
ylabel('$\eta \times 10^{-5}$ Recovered','Interpreter', 'latex')
legend('Injections', 'Unity Slope, Offset $\eta=5.6\times10^{-6}$','Interpreter', 'latex')
set(gca,'FontSize',16);
set(ll,'LineWidth',2);
set(ll,'MarkerSize',16);
set(l,'LineWidth',2);
grid on
xlim([-10 10])
ylim([-10 10])
xticks(inj*1e5)
yticks(inj*1e5)

%% Save plots

if(true)
    fig2=figure(1);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_Injection.pdf','-dpdf','-r1200')
    
end
