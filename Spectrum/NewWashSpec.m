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

if (false)

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
 
end

% Length of days
lenDays = ceil((timFit(end)-timFit(1))/24/3600);

% Calculate complex torque amplitude
torqFit = P.*(C+i*S);

% Fit parameters
fDay = 1/sidDay; % Daily frequency (Hz)
wDay = 2*pi*fDay; % Daily frequency (rad*Hz)
daySamples = floor(sampF/fDay); % Samples in a day
numDaysFit = 2; % Length of cuts (days)
lenMin = 0*numDaysFit; % Minimum length of cut (samples)

% Remove insanely large chi-squareds
Uraw = Uraw(find(Uraw<20));

%% Galaxy Fits

F = linspace(1/sidDay/10,10/sidDay,500);

% Create vectors
A = [];

for freq = F
    
    numDaysFit = 2/freq/sidDay;
    
    amp = [];
    for index = 0:floor(lenDays/numDaysFit)
    
        % Find cut indices
        indexCut = find(floor((timFit-timFit(1))/24/3600/numDaysFit)==index);
    
        % Cut vectors
        cut = timFit(indexCut);
        pol = sign(mean(P(indexCut)));
        u = U(indexCut);
        y = torqFit(indexCut)-mean(torqFit(indexCut));    
            
        if not(isempty(cut))   
            if (length(y)>=lenMin)               
                
                % Design matrix
                x = [cos(2*pi*freq*cut')+i*sin(2*pi*freq*cut')];
    
                if not(isempty(x))
    
                    % Linear least squares fitting to basis functions
                    a = inv(x'*x)*x'*y';
                
                    % Append valid data points
                    if ((a(1)~=0))
                        amp = [amp a(1)];                        
            
                    end
                end
            end
        end
    end
    A = [A abs(mean(amp))];

end

%% Figures


figure(1)
ll = loglog(F,A/m/r,[1/sidDay 1/sidDay], [1e-18 1e-15],'k--');
hold on
text(1.1e-5,2e-18,'Sidereal Frequency','FontSize',16,'Interpreter', 'latex','Rotation',90)
hold off
xlabel('Frequency (Hz)','Interpreter', 'latex')
ylabel('Acceleration (m/s$^2$)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(ll,'LineWidth',2);
grid on

%% Save plots

if(false)
    fig2=figure(1);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_Spec.pdf','-dpdf','-r1200')
end
