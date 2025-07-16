warning('off')

%% Parameters

w0 = 2*pi*6.8944e-4; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
Q = 2e5; % Quality factor
kappa = I*w0^2; % Spring constance (N m/rad)
kb = 1.38064852e-23; % Boltzmann's constant (J/K)
T = 293; % Temperature (K)
thetaCalib = 0.0012523; % Autocollimator calibration (rad/(Diff/Sum))
TTFreq = 0.457120e-3; % Turn table frequency (Hz)

deltaSpeed = 1e-6;

%% Data loading

if (false)
    
    % Run number
    run = ['run6935'];

    % Load vectors form tdms
    inDiff = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Diff");
    inSum = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Sum");
    inTim = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Time");
   
    % Flatten vectors
    inDiff = table2array(inDiff{1});
    inSum = table2array(inSum{1});
    inTim = table2array(inTim{1});
    
end

%% Calibration

% Calculate theta from Diff/Sum
inTheta = thetaCalib*inDiff./inSum;

% Sampling frequency
sampF = 1/(inTim(2)-inTim(1));

% Time indices
startIndex = 1;
endIndex = 2/TTFreq*sampF;

% tim vectors
tim = (startIndex:endIndex)'*sampF;
theta = inTheta(startIndex:endIndex);

%%
% Fit parameters
fFit = TTFreq;
wFit = 2*pi*fFit;

% Design matrix
    x = [cos(wFit*tim) sin(wFit*tim)...
        cos(2*wFit*tim) sin(2*wFit*tim)...
        cos(3*wFit*tim) sin(3*wFit*tim)...
        cos(4*wFit*tim) sin(4*wFit*tim)...
        cos(5*wFit*tim) sin(5*wFit*tim)...
        cos(w0*tim) sin(w0*tim)...
        cos(2*w0*tim) sin(2*w0*tim)...
        cos(3*w0*tim) sin(3*w0*tim)...
        cos(4*w0*tim) sin(4*w0*tim)...
        cos(5*w0*tim) sin(5*w0*tim)...
        ones(length(tim),1)];
    
y = theta;
% Linear least squares fitting to basis functions
a = inv(x'*x)*x'*y;

res = a'*x'-y';

%% Figures

if (true)

    % Angle time series
    figure(1)
    l=plot(tim,theta);
    ylabel('Angle (rad)','Interpreter', 'latex')
    xlabel('Time (s)','Interpreter', 'latex')
    set(gca,'FontSize',16);
    set(l,'LineWidth',1.5);
    grid on

    % Angle time series
    figure(2)
    l=plot(tim,res);
    ylabel('Residual Angle (rad)','Interpreter', 'latex')
    xlabel('Time (s)','Interpreter', 'latex')
    set(gca,'FontSize',16);
    set(l,'LineWidth',1.5);
    grid on    

    figure(3)
    l=plot(theta,res);
    ylabel('Residual Angle (rad)','Interpreter', 'latex')
    xlabel('Theta Angle (rad)','Interpreter', 'latex')
    set(gca,'FontSize',16);
    set(l,'LineWidth',1.5);
    grid on   
   
end

