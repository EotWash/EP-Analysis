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

if (false)
    
    % Run number
    run = ['run6972'];

    % Load vectors form tdms
    inTim = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Time");
    inCycle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="CycleMark");
    inMag = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="BearingPressure");
   
    % Flatten vectors
    inTim = table2array(inTim{1});
    inCycle = table2array(inCycle{1});
    inMag = table2array(inMag{1});
    
    inTim = mod(inTim,31556926)+31556926;

end

if(false)
    rawGal=load('Basis Functions\galVectMin.out');
    galSampF = 1/(rawGal(2,1)-rawGal(1,1))/3600/24;
    timGal=resample(rawGal(:,1),1,ceil(galSampF/sampF));
    inGal=resample(rawGal(:,2),1,ceil(galSampF/sampF));
    outGal=detrend(resample(rawGal(:,3),1,ceil(galSampF/sampF)));
end

%% Calibration


% Sampling frequency
sampF = 1/(inTim(2)-inTim(1));

% Time indices
startIndex = 1;
endIndex = length(inTim);
% startIndex = 8e4;
% endIndex = 8e4;

% Cut vectors
tim = inTim(startIndex:endIndex);
mag = 1e-3*inMag(startIndex:endIndex);
Cycle = inCycle(startIndex:endIndex);

% 1 mHz low pass to remove autocollimator noise
[b,a] = butter(2,2*1e-3/sampF,'low');
mag = filter(b,a,mag);
trim = 2e3;
tim = tim(trim:end);
mag = mag(trim:end);

%%


% Vector creation
C = [];
S = [];

numDaysFit = 2; % Length of cuts (days)
% Length of days
lenDays = ceil((tim(end)-tim(1))/24/3600);

% Cut base fitting
for index = 0
    
    indexCut = find(floor((tim-tim(1))/24/3600/numDaysFit)==index);

    % Cut time vector
    cut = tim(indexCut)';

    % Sync basis function and data 
    galIndex = [];
    for cutGal = cut
        [m,minI] = min(abs(timGal-cutGal/24/3600));
        galIndex = [galIndex minI];
    end
%             
    % Design matrix
    x = [inGal(galIndex)+i*outGal(galIndex)];
    
    % Cut torque vecotor
    y = detrend(mag(indexCut));

    % Linear least squares fitting to basis functions
    a = inv(x'*x)*x'*y;

    % Extract relavent fit parameters
    C = [C real(a(1))];
    S = [S imag(a(1))];        
    
end

mu = mean(C);
un = std(C)/sqrt(length(C));

magAmp = 1e-6;

eta = mu*magAmp/(r*m*aGalaxy);
etaU = un*magAmp/(r*m*aGalaxy);

disp(['Magnetic Amplitude: ' num2str(mu,3) ' G +- ' num2str(un,2) 'G (' num2str(100*un/abs(mu),2) '%)'])

%%
nAv = 1;
[Am, Fm] = asd2(mag, 1/sampF, nAv, 3, @hann);

%% Figures


% Angle time series
figure(1)
l=plot(tim, detrend(mag), cut, real(a'*x'));
ylabel('Magnetic Field (G)','Interpreter', 'latex')
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
l=loglog(Fm,Am);
ylabel('Magnetic Field (G/$\sqrt{Hz}$)','Interpreter', 'latex')
xlabel('Frequency (Hz)','Interpreter', 'latex')
% text(1.5e-4,8e-14,['Magnetic Amplitude: ' num2str(mu,2) ' +- ' num2str(un,2) ' Nm/G'],'FontSize',16,'Interpreter', 'latex')
legend('Data', 'Magnetic', 'TT Frequency','Resonance','Thermal','Interpreter', 'latex')
title('Magnetic Inj','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
% ylim([1e-19 1e-12])
% xlim([1e-5 1e-1])
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
    print(fig2,'EP_MagSpec.pdf','-dpdf','-r1200')

end

