warning('off')

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
aGalaxy =  5e-11; % Acceleration towards dark matter at center of Galaxy (m/s^2)

%% Data loading

if (true)
    
    % Run number
    run = ['run6973'];

    % Load vectors form tdms
    inDiff = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Diff");
    inSum = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Sum");
    inTim = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Time");
    inCycle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="CycleMark");
    inAGI = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="AGI1");
    inAGI2 = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="AGI2");
   
    % Flatten vectors
    inDiff = table2array(inDiff{1});
    inSum = table2array(inSum{1});
    inTim = table2array(inTim{1});
    inCycle = table2array(inCycle{1});
    inAGI = table2array(inAGI{1});
    inAGI2 = table2array(inAGI2{1});
    
end

%% Calibration

% Calculate theta from Diff/Sum
inTheta = inDiff./inSum*thetaCalib;

% Sampling frequency
sampF = 1/(inTim(2)-inTim(1));

% Time indices
% startIndex = 1;
endIndex = length(inTim);
startIndex = 1.45e4;
endIndex = 7.8e4;

% Cut vectors
tim = (startIndex:endIndex)*sampF;
theta = inTheta(startIndex:endIndex);
agi = 1e-5*inAGI(startIndex:endIndex);
agi2 = 1e-5*inAGI2(startIndex:endIndex);
Cycle = inCycle(startIndex:endIndex);

% Torsional response inversion filter
filt = zpk(-2*pi*[pairQ(w0/2/pi,Q)],-2*pi*[pairQ(0.3,1) pairQ(0.3,1)],1);
filt = filt/abs(freqresp(filt,2*pi*1e-4));

% Torque from angle calculation
torqFilt = kappa*lsim(filt, theta, tim);

% Set time zero to first cycle mark
% torqFilt = torqFilt(find(diff(Cycle)>0.5,1):end);
% agi = agi(find(diff(Cycle)>0.5,1):end);
% timFilt = tim(find(diff(Cycle)>0.5,1):end);
timFilt = tim;

% 1 mHz low pass to remove autocollimator noise
[b,a] = butter(2,2*1e-3/sampF,'low');
tqFit = filter(b,a,torqFilt);
agiFit = filter(b,a,agi);
trim = 2e3;
tqFit = tqFit(trim:end);
tFit = timFilt(trim:end);
agiFit = agiFit(trim:end);

% Fit parameters
fFit = TTFreq;
wFit = 2*pi*fFit;
wTilt = 5/4*wFit;
fitSamples =4*floor(sampF/(wTilt/2/pi));

%%
x = [cos(wFit*tFit') sin(wFit*tFit')...
    cos(2*wFit*tFit') sin(2*wFit*tFit')...
    cos(3*wFit*tFit') sin(3*wFit*tFit')...
    cos(4*wFit*tFit') sin(4*wFit*tFit')...
    cos(5*wFit*tFit') sin(5*wFit*tFit')...
    cos(w0*tFit') sin(w0*tFit')...
    cos(2*w0*tFit') sin(2*w0*tFit')...
    cos(3*w0*tFit') sin(3*w0*tFit')...
    cos(4*w0*tFit') sin(4*w0*tFit')...
    cos(5*w0*tFit') sin(5*w0*tFit')...
    ones(length(tFit'),1)];

% Cut torque vecotor
y = tqFit;

% Linear least squares fitting to basis functions
a = inv(x'*x)*x'*y;

% tqFit = tqFit -(a'*x')';

% aa = inv(x'*x)*x'*agiFit;

% agiFit = agiFit -(aa'*x')';
    %%


% Vector creation
C = [];
S = [];

% Cut base fitting
for index = 0:floor(length(tFit)/fitSamples)-1
    
    % Cut time vector
    cut = detrend(tFit(index*fitSamples+1:(index+1)*fitSamples+1))';
    
    % Design matrix
    x = [detrend(agiFit(index*fitSamples+1:(index+1)*fitSamples+1))];
    
    % Cut torque vecotor
    y = tqFit(index*fitSamples+1:(index+1)*fitSamples+1);

    % Linear least squares fitting to basis functions
    a = inv(x'*x)*x'*y;
    
    % Extract relavent fit parameters
    C = [C a(1)];
    S = [S 0];        

end

mu = mean(C);
un = std(C)/sqrt(length(C));

tiltAmp = 4.7e-9;

eta = mu*tiltAmp/(r*m*aGalaxy);
etaU = un*tiltAmp/(r*m*aGalaxy);

disp(['Tilt Coupling: ' num2str(mu,3) ' Nm/rad +- ' num2str(un,2) ' Nm/rad (' num2str(100*un/abs(mu),2) '%)'])
disp(['Tilt Coupling: ' num2str(mu/kappa,3) ' rad/rad +- ' num2str(un/kappa,2) ' rad/rad (' num2str(100*un/abs(mu),2) '%)'])
disp(['eta Tilt: ' num2str(eta,3) ' +- ' num2str(etaU,2)])

%%
nAv = 1;
[A, F] = asd2(theta, 1/sampF, nAv, 3, @hann);
[Am, Fm] = asd2(agiFit, 1/sampF, nAv, 3, @hann);
[Aff, Fff] = asd2(tqFit, 1/sampF, nAv, 3, @hann);

%% Figures


% Angle time series
figure(1)
l=plot(tFit,detrend(tqFit,'constant'), tFit, agiFit*mu);
ylabel('Angle (rad)','Interpreter', 'latex')
xlabel('Time (s)','Interpreter', 'latex')
% ylim([0 1.4e-3])
% xlim([300 640])
%     text(200,1.25e-3, '1 $\mu$Hz Speed Change at 499 s ','FontSize',16,'Interpreter', 'latex')
%     text(200,1.15e-3, txt,'FontSize',16,'Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
grid on

% Torque ASD
figure(4)
l=loglog(Fff,Aff,Fm,Am*abs(mu),[TTFreq TTFreq], [1e-19 1e-14],'--',[w0/2/pi w0/2/pi], [1e-19 1e-14],'--');
ylabel('Torque (N m/$\sqrt{Hz}$)','Interpreter', 'latex')
xlabel('Frequency (Hz)','Interpreter', 'latex')
text(1.5e-4,8e-14,['Tilt Coupling: ' num2str(abs(mu)/kappa*100,2) ' +- ' num2str(un/kappa*100,2) ' %'],'FontSize',16,'Interpreter', 'latex')
legend('Data', 'Coupling $\times$ Tilt', 'TT Frequency','Resonance','Thermal','Interpreter', 'latex')
title('Tilt Injection Foot A 180','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
ylim([1e-19 1e-12])
xlim([1e-5 1e-1])
grid on

figure(5)
l=plot(C*1e12);
ylabel('Coupling Fit (pN m/rad)','Interpreter', 'latex')
xlabel('Cut Number','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(l,'LineWidth',1.5);
% xlim([-0.25 0.25])
% ylim([-0.25 0.25])
grid on

    %%
if (true)
    fig2=figure(4);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_TiltInjFootA180.pdf','-dpdf','-r1200')

end

