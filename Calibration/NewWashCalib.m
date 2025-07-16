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

deltaSpeed = 1e-6;

%% Data loading

if (true)
    
    % Run number
    run = ['run6932'];

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
inTheta = inDiff./inSum*thetaCalib;

% Sampling frequency
sampF = 1/(inTim(2)-inTim(1));

% Time indices
startIndex = 1;
endIndex = length(inTim);

% Cut vectors
tim = (startIndex:endIndex)*sampF;
theta = inTheta(startIndex:endIndex);

% Time indices
startIndex = 500;
endIndex = 550;

timM = (startIndex:endIndex)*sampF;
thetaM = inTheta(startIndex:endIndex)/thetaCalib;

dAng=mean(diff(thetaM')./diff(timM));
dAngU=std(diff(thetaM')./diff(timM));
thetaCalibM = deltaSpeed*2*pi/dAng;
thetaCalibU = deltaSpeed*2*pi/dAng^2*dAngU/sqrt(length(thetaM));

txt = ['Calibration: ' num2str(thetaCalibM,3) ' +- ' num2str(thetaCalibU,2) ' (' num2str(100*thetaCalibU/thetaCalib,2) '%)'];
disp(txt)

%% Figures


% Angle time series
figure(1)
l=plot(tim,theta, tim , 2*pi*deltaSpeed*tim+(2.6e-4-2*pi*deltaSpeed*499));
legend('Data','Expected from 1 $\mu$Hz Speed Change','Interpreter', 'latex','Location','northwest')
ylabel('Angle (rad)','Interpreter', 'latex')
xlabel('Time (s)','Interpreter', 'latex')
ylim([0 1.4e-3])
xlim([300 640])
%     text(200,1.25e-3, '1 $\mu$Hz Speed Change at 499 s ','FontSize',16,'Interpreter', 'latex')
%     text(200,1.15e-3, txt,'FontSize',16,'Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
grid on
    %%
if (true)
    fig2=figure(1);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_Calibration.pdf','-dpdf','-r1200')

end

