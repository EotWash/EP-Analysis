warning('off')

%% Data loading

if (false)

    run = 'run6875';

    inTTAngle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Angle");
    inDiff = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Diff");
    inSum = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Sum");
    inTim = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="Time");
    inCycle = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="CycleMark");
    inPress = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="IonPump");
    inAGI = tdmsread(['G:\Shared drives\Eot-Wash\NewWash\Data\' run '.tdms'], ChannelGroup="raw_data", ChannelNames="AGI1");

    inTTAngle = table2array(inTTAngle{1});
    inDiff = table2array(inDiff{1});
    inSum = table2array(inSum{1});
    inTim = table2array(inTim{1});
    inCycle = table2array(inCycle{1});
    inPress = table2array(inPress{1});
    inAGI = table2array(inAGI{1});

    inDEP = load('DarkEPData.csv');
    fDEP = inDEP(:,1);
    aDEP = inDEP(:,2);
end

inTheta = inDiff./inSum;

sampF = 1/(inTim(2)-inTim(1));

% startIndex = 1.04e5;
% startIndex = 6.5e4;
startIndex = 10;
endIndex = length(inTim);
% endIndex = 4e;
% endIndex = 1.2e5;

%% Calibration

w0 = 2*pi*6.8944e-4;
I = 3.78e-5;
Q = 2e5;
kappa = I*w0^2;
% kappa = 7.129e-10; % Theoretical 
kb = 1.38064852e-23;
T = 293;
thetaCalib = 3/300/8;
m = 20e-3;
r = 2e-2;

Msun = 1.9891e30; % Mass of sun (kg)
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
Rsun = 149.6e9; % Radius from earth to sun (m)

aSun = G*Msun/Rsun^2; % Acceleration towards sun (m/s^2)
aEarth = 1.68e-2;

% tim = inTim(startIndex:endIndex);
tim = (startIndex:endIndex)*sampF;
theta = thetaCalib*inTheta(startIndex:endIndex);
TTAngle = 60*inTTAngle(startIndex:endIndex);
Cycle = inCycle(startIndex:endIndex);
Press = inPress(startIndex:endIndex);

% TTFreq = mean(diff(TTAngle/360))*sampF;
TTFreq = 0.457120e-3;

aDM = 1e-25/9.73e-19 % Torque to g_dm conversion for Be-Al


%% Torsional filter

filt = zpk(-2*pi*[pairQ(w0/2/pi,Q)],-2*pi*[pairQ(0.3,1) pairQ(0.3,1)],1);
filt = filt/abs(freqresp(filt,2*pi*1e-4));

torqFilt = kappa*lsim(filt, theta, tim);

torqFilt = torqFilt(find(diff(Cycle)>0.5,1):end);
timFilt = tim(find(diff(Cycle)>0.5,1):end);
angFilt = theta(find(diff(Cycle)>0.5,1):end);
timFilt = timFilt - timFilt(1);

%% Turn table correction fits
fFit = TTFreq;
wFit = 2*pi*fFit;

x = [cos(wFit*timFilt') sin(wFit*timFilt')...
    cos(2*wFit*timFilt') sin(2*wFit*timFilt')...
    cos(3*wFit*timFilt') sin(3*wFit*timFilt')...
    cos(4*wFit*timFilt') sin(4*wFit*timFilt')...
    cos(5*wFit*timFilt') sin(5*wFit*timFilt')...
    cos(w0*timFilt') sin(w0*timFilt')...
    cos(2*w0*timFilt') sin(2*w0*timFilt')...
    cos(3*w0*timFilt') sin(3*w0*timFilt')...
    cos(4*w0*timFilt') sin(4*w0*timFilt')...
    cos(5*w0*timFilt') sin(5*w0*timFilt')...
    timFilt' timFilt'.^2 timFilt'.^3 ...
    ones(length(timFilt'),1)];

y = detrend(angFilt);
a = inv(x'*x)*x'*y;
res = y'-a'*x';

[Ar, Fr] = asd2(res, 1/sampF, 1, 3, @hann);

figure(111)
plot(timFilt, res)

figure(112)
loglog(Fr, Ar)

disp(['Cosine 1omega: ' num2str(a(1))])
disp(['Sine 1omega: ' num2str(a(2))])
disp(['Cosine 2omega: ' num2str(a(3))])
disp(['Sine 2omega: ' num2str(a(4))])
disp(['Cosine 3omega: ' num2str(a(5))])
disp(['Sine 3omega: ' num2str(a(6))])
disp(['Cosine 4omega: ' num2str(a(7))])
disp(['Sine 4omega: ' num2str(a(8))])
%% ASD Calculation

nAv = 1;
[A, F] = asd2(theta, 1/sampF, nAv, 3, @hann);
[Af, Ff] = asd2(torqFilt, 1/sampF, nAv, 3, @hann);

%% Torque Calculation
w = 2*pi*F;

R = 1./(1-w.^2/w0.^2-i/Q)/kappa; %% Torq to Angle Response

torq = abs(A./R');

%% Themal Limit

thermT = abs(sqrt(4*kb*T*(kappa/Q).*(1./w)));
thermA = abs(R.*thermT);

%% Fits
fFit = TTFreq;
wFit = 2*pi*fFit;
fitSamples = 2*floor(sampF/fFit);

[b,a] = butter(3,2*1e-3/sampF,'low');
tqFit = filter(b,a,torqFilt);
tqFit = tqFit(2e3:end);
tFit = timFilt(2e3:end);
ttFit = TTAngle(2e3:end);

[Aff, Fff] = asd2(tqFit, 1/sampF, nAv, 3, @hann);

C = [];
S = [];
U = [];
CR = [];
SR = [];
timFit = [];

for index = 0:floor(length(tFit)/fitSamples)-1
    
    cut = tFit(index*fitSamples+1:(index+1)*fitSamples+1)';

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
    
    y = tqFit(index*fitSamples+1:(index+1)*fitSamples+1);

    a = inv(x'*x)*x'*y;
    
    C = [C a(1)];
    S = [S a(2)];
    CR = [CR a(11)];
    SR = [SR a(12)];

    U = [U std(a'*x'-y')];

    timFit = [timFit mean(timFilt(index*fitSamples+1:(index+1)*fitSamples+1))];

end

%% Uncertainty cuts

unCut = find(U<mean(U));
C = C(unCut);
S = S(unCut);
U = U(unCut);
timFit = timFit(unCut);
CR = CR(unCut);
SR = SR(unCut);

%% Q Calculations

AR = sqrt(CR.^2+SR.^2);

x = [timFit; ones(1,length(timFit))]';
        
y = log(AR');

a = inv(x'*x)*x'*y;

QMeas = w0/2/a(1);

disp(['Q Measured: ' num2str(QMeas)])
disp([' '])

%%

torqFit = C+i*S;

amp = sqrt(mean(C)^2+mean(S)^2);
unc = 1/amp*sqrt(mean(C)^2*std(C)^2+mean(S)^2*std(S)^2)/sqrt(length(C));

etaEarth = amp/(r*m*aEarth);
etaEarthUnc = unc/(r*m*aEarth);

disp(['Cosine Torque: ' num2str(mean(C)*1e15) ' fN m +- ' num2str(std(C)/sqrt(length(C))*1e15) ' fN m'])
disp(['Sine Torque: ' num2str(mean(S)*1e15) ' fN m +- ' num2str(std(S)/sqrt(length(S))*1e15) ' fN m'])
disp(['Amp Torque: ' num2str(amp*1e15) ' fN m +- ' num2str(unc*1e15) ' fN m'])
disp([' '])
disp(['Eta Earth: ' num2str(etaEarth) ' +- ' num2str(etaEarthUnc)])
disp([' '])

%% Thermal Circle Calculations

thermPhi = linspace(0,2*pi,100);
thermAmp = abs(sqrt(4*kb*T*(kappa/Q).*(1./(2*pi*TTFreq))))*sqrt((2*pi*TTFreq));
thermCirc = thermAmp*(cos(thermPhi)+i*sin(thermPhi))+mean(torqFit);

%% Dark Matter Calcuations

amp1w = sqrt(C.^2+S.^2);

sampFDM = sampF/fitSamples;
[ADM, FDM] = asd2(amp1w, 1/sampFDM, nAv, 3, @hann);
ADM = aDM*ADM/sqrt(length(C)/sampFDM);
ADMT = aDM*Aff/sqrt(timFilt(end));

fDay = 1/24/3600;
wDay = 2*pi*fDay;

dmFreq = linspace(1/timFilt(end),0.95*sampFDM/2,floor(sampFDM*timFilt(end)));

ampDM = [];
uncDM = [];

for indexDM = 1:length(dmFreq)
    
    fFitDM = dmFreq(indexDM);
    wFit = 2*pi*fFitDM;
    fitSamplesDM = 2*floor(sampFDM/fFitDM);
    
    CDM = [];
    SDM = [];
    
    for index = 0:floor(length(timFit)/fitSamplesDM)-2
        
        cut = timFit(index*fitSamplesDM+1:(index+1)*fitSamplesDM+1)';
  
%         x = [cos(wFit*cut).*cos(wDay*cut) cos(wFit*cut).*sin(wDay*cut) ...
%             sin(wFit*cut).*cos(wDay*cut) sin(wFit*cut).*sin(wDay*cut)...        
%             ones(length(cut),1)];
        x = [cos(wFit*cut) sin(wFit*cut)...        
            ones(length(cut),1)];
        
        y = amp1w(index*fitSamplesDM+1:(index+1)*fitSamplesDM+1)';
    
        a = inv(x'*x)*x'*y;
        
%         CDM = [CDM sqrt(a(1)^2+a(2)^2)];
%         SDM = [SDM sqrt(a(3)^2+a(4)^2)];
        CDM = [CDM a(1)];
        SDM = [SDM a(2)];
    
    end
        
    ampDM = [ampDM; sqrt(mean(CDM)^2+mean(SDM)^2)];
    uncDM = [uncDM; 1/amp*sqrt(mean(CDM)^2*std(CDM)^2+mean(SDM)^2*std(SDM)^2)/sqrt(length(CDM))];

end

ampDM = aDM*ampDM;
uncDM = aDM*uncDM;


%% Figures

figure(1)
l=plot(tim,theta);
ylabel('Angle (rad)','Interpreter', 'latex')
xlabel('Time (s)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
grid on

mF = logspace(-5,-2);
mA = mF*0+1e-25;

figure(10)
l=loglog(FDM, ADM, mF, mA, fDEP, aDEP, Fff, ADMT);
ylabel('$g_{DM}$','Interpreter', 'latex')
xlabel('Frequency (Hz)','Interpreter', 'latex')
legend('Data', 'MICROSCOPE','DarkEP','Raw Torq Limits','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
ylim([1e-27 1e-23])
xlim([1e-5 1e-2])
grid on

figure(12)
l=loglog(dmFreq, ampDM, mF, mA, '--', fDEP, aDEP, '--');
ylabel('$g_{DM}$','Interpreter', 'latex')
xlabel('Frequency (Hz)','Interpreter', 'latex')
legend('Data', 'MICROSCOPE','DarkEP','Raw Torq Limits','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
ylim([1e-27 1e-23])
xlim([0.9e-5 1e-3])
grid on

if(false)
    fig2=figure(12);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'ULDM_Limits.pdf','-dpdf','-r1200')
    
   
end
