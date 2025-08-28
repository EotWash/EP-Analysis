warning('off')

%% Parameters

w0 = 2*pi/1461; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
Q = 1.1e5; % Quality factor
kappa = I*w0^2; % Spring constance (N m/rad)
kb = 1.38064852e-23; % Boltzmann's constant (J/K)
T = 293; % Temperature (K)
thetaCalib = 3/300/8; % Autocollimator calibration (rad/(Diff/Sum))
TTFreq = 0.4568e-3; % Turn table frequency (Hz)

%% Data loading

if (true)
    
    % Run number
    run = ['run6964'];

    % Load vectors form tdms
    inTTAngle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Angle");
    inDiff = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Diff");
    inSum = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Sum");
    inTim = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Time");
    inCycle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="CycleMark");
   
    % Flatten vectors
    inTTAngle = table2array(inTTAngle{1});
    inDiff = table2array(inDiff{1});
    inSum = table2array(inSum{1});
    inTim = table2array(inTim{1});
    inCycle = table2array(inCycle{1});
    
end

%% Calibration

% Calculate theta from Diff/Sum
inTheta = thetaCalib*inDiff./inSum;

% Sampling frequency
sampF = 1/(inTim(2)-inTim(1));

% Time indices
startIndex = 1e1;
endIndex =length(inTTAngle);

% Cut vectors
tim = (startIndex:endIndex)*sampF;
theta = inTheta(startIndex:endIndex);
TTAngle = 60*inTTAngle(startIndex:endIndex);
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
[A, F] = asd2(theta, 1/sampF, nAv, 3, @hann);
[Af, Ff] = asd2(torqFilt, 1/sampF, nAv, 3, @hann);
[Aff, Fff] = asd2(tqFit, 1/sampF, nAv, 3, @hann);

%% Themal Limit

w = 2*pi*F;

R = 1./(1-w.^2/w0.^2-i/Q)/kappa; %% Torq to Angle Response

thermT = abs(sqrt(4*kb*T*(kappa/Q).*(1./w)));
thermA = abs(R.*thermT);

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
    
    % Design matrix
    x = [cos(wFit*cut) sin(wFit*cut)...
        cos(2*wFit*cut) sin(2*wFit*cut)...
        cos(3*wFit*cut) sin(3*wFit*cut)...
        cos(4*wFit*cut) sin(4*wFit*cut)...
        cos(5*wFit*cut) sin(5*wFit*cut)...
        cos(w0*cut) sin(w0*cut)...
        cos(2*w0*cut) sin(2*w0*cut)...
        cos(3*w0*cut) sin(3*w0*cut)...
        cos(4*w0*cut) sin(4*w0*cut)...
        cos(5*w0*cut) sin(5*w0*cut)...
        ones(length(cut),1)];
    
    % Cut torque vecotor
    y = tqFit(index*fitSamples+1:(index+1)*fitSamples+1);

    % Linear least squares fitting to basis functions
    a = inv(x'*x)*x'*y;
    
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

if (true)

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

