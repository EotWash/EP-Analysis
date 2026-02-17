% warning('off')

%% Parameters

w0 = 2*pi/1461; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
Q = 1.1e5; % Quality factor
kappa = I*w0^2; % Spring constance (N m/rad)
kb = 1.38064852e-23; % Boltzmann's constant (J/K)
T = 293; % Temperature (K)
thetaCalib = 3/300/8; % Autocollimator calibration (rad/(Diff/Sum))

%% Data loading
runs = ["run6891" "run6893" ...
   "run6895" "run6896" "run6897" "run6900" "run6903" "run6904" "run6905" "run6923" "run6925" "run6926" "run6927" "run6930"...
    "run6931" "run6936" "run6939" "run6949" "run6950" "run6954" "run6955" "run6956" "run6958" "run6962" "run6964"]; % Data pre-break
runs=[runs,"run6979", "run6981", "run6982" "run6984", "run6985", "run6986", "run6987", "run6988"]; %post break data

for j=1:length(runs) %Looping for all runs

if (true)
    
    % Run number
    run = runs{j} %set which run

    % Load vectors form tdms

    %Michael's path
    % inTTAngle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Angle");
    % inDiff = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Diff");
    % inSum = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Sum");
    % inTim = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Time");
    % inCycle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="CycleMark");
    % setFreq = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="SetFrequency"); %frequency which TT set to

    %Shoshana's Path
    inTTAngle = tdmsread(['H:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Angle");
    inDiff = tdmsread(['H:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Diff");
    inSum = tdmsread(['H:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Sum");
    inTim = tdmsread(['H:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Time");
    inCycle = tdmsread(['H:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="CycleMark");
    setFreq = tdmsread(['H:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="SetFrequency"); %frequency which TT set to

    % Flatten vectors
    inTTAngle = table2array(inTTAngle{1});
    inDiff = table2array(inDiff{1});
    inSum = table2array(inSum{1});
    inTim = table2array(inTim{1});
    inCycle = table2array(inCycle{1});
    setFreq=table2array(setFreq{1});
    TTFreq=mode(round(setFreq,4))*1e-3; %Finds what the TTFreq is set to
end

%% Calibration

% Calculate theta from Diff/Sum
inTheta = thetaCalib*inDiff./inSum;

% Sampling frequency
sampF = 1/mode(round(diff(inTim),3)); %Getting sample frequency

% Time indices
startIndex = 1e1;
endIndex =length(inTTAngle);

% Cut vectors
tim = (startIndex:endIndex)*sampF;
theta = inTheta(startIndex:endIndex);
TTAngle = inTTAngle(startIndex:endIndex);
Cycle = inCycle(startIndex:endIndex);

%% Torsional filter

% Torsional response inversion filter
filt = zpk(-2*pi*[pairQ(w0/2/pi,Q)],-2*pi*[pairQ(0.3,1) pairQ(0.3,1)],1);
filt = filt/abs(freqresp(filt,2*pi*1e-4));

% Torque from angle calculation
torqFilt = kappa*lsim(filt, theta, tim);

% Set time zero to first cycle mark
torqFilt = torqFilt(find(diff(Cycle)>0.5,1):end);
timFilt = tim(find(diff(Cycle)>0.5,1):end);
angFilt = theta(find(diff(Cycle)>0.5,1):end);

% 1 mHz low pass to remove autocollimator noise
[b,a] = butter(3,2*1e-3/sampF,'low');
tqFit = filter(b,a,torqFilt);
tqFit = tqFit(2e3:end);
tFit = timFilt(2e3:end);
ttFit = TTAngle(2e3:end);

%% ASD Calculation

nAv = 1;
% [A, F] = asd2(theta, 1/sampF, nAv, 3, @hann);
% [Af, Ff] = asd2(torqFilt, 1/sampF, nAv, 3, @hann);
% [Aff, Fff] = asd2(tqFit, 1/sampF, nAv, 3, @hann);
% 
% %% Themal Limit
% 
% w = 2*pi*F;
% 
% R = 1./(1-w.^2/w0.^2-1i/Q)/kappa; %% Torq to Angle Response
% 
% thermT = abs(sqrt(4*kb*T*(kappa/Q).*(1./w)));
% thermA = abs(R.*thermT);

%% Fits

% Fit parameters
fFit = TTFreq;
wFit = 2*pi*fFit;
fitSamples = 2*floor(sampF/fFit);

% Vector creation
C = [];
S = [];
U = [];
CR = [];
SR = [];
timFit = [];


% Cut base fitting
for index = 0:floor(length(tFit)/fitSamples)-1
    
    % Cut time vector
    cut = tFit(index*fitSamples+1:(index+1)*fitSamples+1)';
    fitAngle=ttFit(index*fitSamples+1:(index+1)*fitSamples+1); %Using angle data to fit
    % Design matrix
    x = [cos(fitAngle) sin(fitAngle)...
        cos(2*fitAngle) sin(2*fitAngle)...
        cos(3*fitAngle) sin(3*fitAngle)...
        cos(4*fitAngle) sin(4*fitAngle)...
        cos(5*fitAngle) sin(5*fitAngle)...
        cos(w0*cut) sin(w0*cut)...
        cos(2*w0*cut) sin(2*w0*cut)...
        cos(3*w0*cut) sin(3*w0*cut)...
        cos(4*w0*cut) sin(4*w0*cut)...
        cos(5*w0*cut) sin(5*w0*cut)...
        ones(length(cut),1)];

    % Cut torque vecotor
    y = tqFit(index*fitSamples+1:(index+1)*fitSamples+1);

    % Linear least squares fitting to basis functions
    a = (x'*x)\x'*y;
    
    % Extract relavent fit parameters
    C = [C a(1)];
    S = [S a(2)];
    CR = [CR a(11)];
    SR = [SR a(12)];
    
    % Calculate misfit
    U = [U std(a'*x'-y')];
    
    % Time stamp
    timFit = [timFit mean(cut)];

end

%% Save fit amplitudes

if (true)
    out = [timFit+inTim(1); C; S; U];
    save([run 'Fits.mat'],'out')
end

% Make complex torque amplitude
torqFit = C+i*S;

%% Thermal Circle Calculations

thermPhi = linspace(0,2*pi,100);
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq));
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

%% Figures
end
if (false)

    % Angle time series
    figure(1)
    l=plot(tim,inDiff(startIndex:endIndex));
    ylabel('Angle (rad)','Interpreter', 'latex')
    xlabel('Time (s)','Interpreter', 'latex')
    set(gca,'FontSize',16);
    set(l,'LineWidth',1.5);
    grid on
    
    % Torque time series
    figure(2)
    l=plot(tFit, tqFit);
    ylabel('Torque (N m)','Interpreter', 'latex')
    xlabel('Time (days)','Interpreter', 'latex')
    set(gca,'FontSize',16);
    set(l,'LineWidth',1.5);
    grid on
    
    % Angle ASD
    figure(3)
    l=loglog(F,A,[TTFreq TTFreq], [1e-8 1e-4],'--',[w0/2/pi w0/2/pi], [1e-8 1e-4],'--',F,thermA);
    ylabel('Angle (rad/$\sqrt{Hz}$)','Interpreter', 'latex')
    xlabel('Frequency (Hz)','Interpreter', 'latex')
    legend('Data', 'TT Frequency' ,'Resonance','Thermal','Interpreter', 'latex')
    set(gca,'FontSize',16);
    set(l,'LineWidth',1.5);
    ylim([1e-12 1e-2])
    xlim([1e-5 1e0])
    grid on
    
    % Torque ASD
    figure(4)
    l=loglog(Fff,Aff,[TTFreq TTFreq], [1e-19 1e-14],'--',[w0/2/pi w0/2/pi], [1e-19 1e-14],'--',F,thermT);
    ylabel('Torque (N m/$\sqrt{Hz}$)','Interpreter', 'latex')
    xlabel('Frequency (Hz)','Interpreter', 'latex')
    legend('Data', 'TT Frequency','Resonance','Thermal','Interpreter', 'latex')
    set(gca,'FontSize',16);
    set(l,'LineWidth',1.5);
    ylim([1e-17 1e-12])
    xlim([1e-5 2e-3])
    grid on    
       
    % Fit quadrature plot
    figure(5)
    l=plot(real(torqFit)*1e15,imag(torqFit)*1e15,'.',real(thermCirc)*1e15,imag(thermCirc)*1e15);
    ylabel('Sine Torque (fN m)','Interpreter', 'latex')
    xlabel('Cosine Torque (fN m)','Interpreter', 'latex')
    legend('Data','Thermal Noise','Interpreter', 'latex')
    set(gca,'FontSize',16);
    set(l,'MarkerSize',16);
    set(l,'LineWidth',1.5);
    grid on
    
    % Fit amplitude time series
    figure(6)
    l=errorbar(timFit/3600/24, abs(torqFit)*1e15, U*1e15,'.');
    ylabel('Torque Amplitude (fN m)','Interpreter', 'latex')
    xlabel('Time (days)','Interpreter', 'latex')
    set(gca,'FontSize',16);
    set(l,'MarkerSize',16);
    grid on
    
end

