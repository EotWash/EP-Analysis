
warning('off')

loadData = true;

%% Parameters

w0 = 2*pi*6.8944e-4; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
Q = 2e5; % Quality factor
kappa = I*w0^2; % Spring constance (N m/rad)
kb = 1.38064852e-23; % Boltzmann's constant (J/K)
T = 293; % Temperature (K)
thetaCalib = 3/300/8; % Autocollimator calibration (rad/(Diff/Sum))
TTFreq = 0.457120e-3; % Turn table frequency (Hz)

m = 38.72e-3/2; % Mass (kg)
r = 3.77e-2/2; % Lever-arm (m)

%% Expected signal
m4 = 1.25;
R = 28.6e-2;

q = 9.35e-9;
Q = 3/4*sqrt(35/2/pi)*m4*R^-5;
G = 6.67430e-11;

torqExp = 4*pi*G*4/9*q*Q;

%% No Q44 Inj


if (loadData)
    
    % Run number
    run = ['run6931'];

    % Load vectors form tdms   
    inDiff = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Diff");
    inSum = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Sum");
    inCycle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="CycleMark");
   
    % Flatten vectors
    inDiff1 = table2array(inDiff{1});
    inSum1 = table2array(inSum{1});
    inCycle1 = table2array(inCycle{1});
    
end

% Calculate theta from Diff/Sum
inTheta = thetaCalib*inDiff1./inSum1;

% Sampling frequency
sampF = 1;

% Time indices
startIndex = 1;
endIndex =length(inDiff1);

% Cut vectors
tim = (startIndex:endIndex);
theta = inTheta(startIndex:endIndex);
Cycle = inCycle1(startIndex:endIndex);

% Torsional response inversion filter
filt = zpk(-2*pi*[pairQ(w0/2/pi,Q)],-2*pi*[pairQ(0.3,1) pairQ(0.3,1)],1);
filt = filt/abs(freqresp(filt,2*pi*1e-4));

% Torque from angle calculation
torqFilt = kappa*lsim(filt, theta, tim);

% Set time zero to first cycle mark
torqFilt = torqFilt(find(diff(Cycle)>0.5,1):end);
timFilt = tim(find(diff(Cycle)>0.5,1):end);

% 1 mHz low pass to remove autocollimator noise
[b,a] = butter(3,2*1e-2/sampF,'low');
tqFit = filter(b,a,torqFilt);
tqFit = tqFit(2e3:end);
tFit = timFilt(2e3:end);

% Fit parameters
fFit = TTFreq;
wFit = 2*pi*fFit;
fitSamples = 2*floor(sampF/fFit);

% Vector creation
Cn = [];
Sn = [];


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
    Cn = [Cn a(7)];
    Sn = [Sn a(8)];        

end

%% With Q44 Inj

if (loadData)
    
    % Run number
    run = ['run6936'];

    % Load vectors form tdms   
    inDiff = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Diff");
    inSum = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Sum");
    inCycle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="CycleMark");
   
    % Flatten vectors
    inDiff2 = table2array(inDiff{1});
    inSum2 = table2array(inSum{1});
    inCycle2 = table2array(inCycle{1});
    
end

% Calculate theta from Diff/Sum
inTheta = thetaCalib*inDiff2./inSum2;

% Time indices
startIndex = 1;
endIndex =length(inDiff2);

% Cut vectors
tim = (startIndex:endIndex);
theta = inTheta(startIndex:endIndex);
Cycle = inCycle2(startIndex:endIndex);


% Torque from angle calculation
torqFilt = kappa*lsim(filt, theta, tim);

% Set time zero to first cycle mark
torqFilt = torqFilt(find(diff(Cycle)>0.5,1):end);
timFilt = tim(find(diff(Cycle)>0.5,1):end);

[b,a] = butter(3,2*1e-2/sampF,'low');
tqFit = filter(b,a,torqFilt);
tqFit = tqFit(2e3:end);
tFit = timFilt(2e3:end);

% Vector creation
C = [];
S = [];

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
    C = [C a(7)];
    S = [S a(8)];        

end

%% With Q44 minus two masses Inj

if (loadData)
    
    % Run number
    run = ['run6939'];

    % Load vectors form tdms   
    inDiff = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Diff");
    inSum = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Sum");
    inCycle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="CycleMark");
   
    % Flatten vectors
    inDiff3 = table2array(inDiff{1});
    inSum3 = table2array(inSum{1});
    inCycle3 = table2array(inCycle{1});
    
end

% Calculate theta from Diff/Sum
inTheta = thetaCalib*inDiff3./inSum3;

% Time indices
startIndex = 1;
endIndex =length(inDiff3);

% Cut vectors
tim = (startIndex:endIndex);
theta = inTheta(startIndex:endIndex);
Cycle = inCycle3(startIndex:endIndex);


% Torque from angle calculation
torqFilt = kappa*lsim(filt, theta, tim);

% Set time zero to first cycle mark
torqFilt = torqFilt(find(diff(Cycle)>0.5,1):end);
timFilt = tim(find(diff(Cycle)>0.5,1):end);

[b,a] = butter(3,2*1e-2/sampF,'low');
tqFit = filter(b,a,torqFilt);
tqFit = tqFit(2e3:end);
tFit = timFilt(2e3:end);

% Vector creation
C2 = [];
S2 = [];

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
    C2 = [C2 a(7)];
    S2 = [S2 a(8)];        

end
%% 

torqMeas = abs((mean(C)-mean(Cn))+i*(mean(S)-mean(Sn)))
per = (torqMeas/torqExp-1)*100

%% Figures


phi = linspace(0,2*pi);
figure(5)
l=plot(Cn*1e12/m/r,Sn*1e12/m/r,'.',C*1e12/m/r,S*1e12/m/r,'.',...
    (torqExp*cos(phi)+mean(Cn))*1e12/m/r,(torqExp*sin(phi)+mean(Sn))*1e12/m/r,'k');
hold on
x = ([0.9*torqExp*cos(phi), 1.1*torqExp*cos(fliplr(phi))]+mean(Cn))*1e12/m/r;
y = ([0.9*torqExp*sin(phi), 1.1*torqExp*sin(fliplr(phi))]+mean(Sn))*1e12/m/r;
fill(x,y,'k','LineStyle','none');
alpha(0.1)
hold off
ylabel('Sine Acceleration (pm/s$^2$)','Interpreter', 'latex')
xlabel('Cosine Acceleration (pm/s$^2$)','Interpreter', 'latex')
legend('No $Q_{4,4}$','With $Q_{4,4}$','Expected Change','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(l,'LineWidth',1.5);
% ylim([-8 4])
% xlim([-4 8])
grid on
%%
if (false)

    fig2=figure(5);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_Q44Inj.pdf','-dpdf','-r1200')

end

