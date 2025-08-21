
galDec = -29-28/3600;

seaLat = 47+36*60+14/3600;
seaLong = 122+19/60+48/3600;

eaR = 6371000;

fDay = 1/24/3600;
wDay = 2*pi*fDay;
sampF = 0.457120e-3/2;

if false
    % Loading in galaxy basis funtions outputted from galVect.py
    rawGal=load('galVectMin.out');
end
galSampF = 1/(rawGal(2,1)-rawGal(1,1))/3600/24;
timGal=decimate(rawGal(:,1),floor(galSampF/sampF))';
inGal=detrend(decimate(rawGal(:,2),floor(galSampF/sampF)));
outGal=detrend(decimate(rawGal(:,3),floor(galSampF/sampF)));

geoGal = 1e10*[sin(galDec*pi/180) 0 cos(galDec*pi/180)];
trans = geoGal + [0 0 eaR];
labGal = trans*[cos((90-seaLat)*pi/180) 0 sin((90-seaLat)*pi/180);
    0 1 0;
    -sin((90-seaLat)*pi/180) 0 cos((90-seaLat)*pi/180)];


inGalAn = [];
for tim = timGal 
    labGalAn = labGal*[cos((2*pi*tim-seaLong)*pi/180) sin((2*pi*tim-seaLong)*pi/180) 0;
    -sin((2*pi*tim-seaLong)*pi/180) cos((2*pi*tim-seaLong)*pi/180) 0;
    0 0 1];
    inGalAn = [inGalAn; labGalAn(2)/1e10];
end

figure(1)
l=plot(timGal,inGal, timGal, inGalAn);%, timGalD, inGalD);
ylabel('Basis function','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
legend('Raw','Decimated')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
grid on
