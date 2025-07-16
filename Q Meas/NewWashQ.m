warning('off')

%% Parameters

w0 = 2*pi*6.8567e-4; % Resonant frequency (rad*Hz)
I = 3.78e-5; % Moment of inertia (kg-m^2)
Q = 2e5; % Quality factor
kappa = I*w0^2; % Spring constance (N m/rad)
kb = 1.38064852e-23; % Boltzmann's constant (J/K)
T = 293; % Temperature (K)
thetaCalib = 0.0012523; % Autocollimator calibration (rad/(Diff/Sum))
TTFreq = 0.457120e-3; % Turn table frequency (Hz)

deltaSpeed = 1e-6;

%% Data loading

if (true)
    
    % Run number
    run = ['run6938'];

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
sampF = 1;

% Time indices
startIndex = 1;
endIndex = length(inTheta);

% tim vectors
tim = (startIndex:endIndex)'*sampF;
theta = detrend(inTheta(startIndex:endIndex));

%% Q Measurement

firstPeak = 1414;
period = floor(1/6.8567e-4)+3;
window = 20;
numPeaks = ceil((tim(end)-firstPeak)/period);


pks = [];
pkt = [];

sumThres = 0.6;
figure(3)
hold off
for index = 0:numPeaks-1
    cut = theta(firstPeak+index*period:firstPeak+index*period+window);
    cutSum = inSum(firstPeak+index*period:firstPeak+index*period+window);
    cut = cut(find(cutSum>sumThres));

    [M,I] = min(cut);
    pks = [pks M];
    pkt = [pkt tim(I)+firstPeak+index*period];

    figure(3)   
    plot(cut)
    hold on
end

theta0 = 90*pi/180;
pks = (theta0 - pks)/theta0;

a = mean(diff(pks)./diff(pkt));
aU = std(diff(pks)./diff(pkt))/sqrt(length(pkt));
b = mean(pks)-a*mean(pkt);

Q = w0/2/abs(a);
QU = w0/2/a^2*aU/sqrt(length(pks));

Q2 = pi*90/abs(pks(end)-pks(1))/length(pks)


% dTheta = mean(diff(pks));



% a = mean(diff(pks)/pks(1));
% aU = std(diff(pks)/pks(1))/sqrt(length(pkt));
% b = mean(pks)-a*mean(pkt);

% Q = pi/abs(a)
% QU = pi/abs(a)^2*aU/sqrt(length(pks));

txt = ['Q: ' num2str(Q,3) ' +- ' num2str(QU,2) ' (' num2str(100*QU/Q,2) '%)'];
disp(txt)
%% Figures


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
l=plot(pkt, pks, '.', pkt ,a*pkt+b);
ylabel('Angle (rad)','Interpreter', 'latex')
xlabel('Time (s)','Interpreter', 'latex')
legend('Data','Fit','Interpreter', 'latex')
text(8000,0.9998, txt,'FontSize',16,'Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
set(l,'MarkerSize',16);
grid on


if (false)
   fig2=figure(2);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_QMeas.pdf','-dpdf','-r1200')

end

