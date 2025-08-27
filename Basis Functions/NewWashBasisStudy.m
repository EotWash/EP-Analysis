
galDec = -(29+28/3600);
galRA = 360/24*(17+45/60+40.0409/3600);

seaLat = 47+36/60+36.94/3600;
seaLong = 122+18/60+13.24/3600;

eaR = 6371000;
galR = 2.5e20;

fDay = 1/24/3600;
wDay = 2*pi*fDay;
sampF = 0.457120e-3/2;

if true
    % Loading in galaxy basis funtions outputted from galVect.py
    rawGal=load('galVectMinTest.out');
end
galSampF = 1/(rawGal(2,1)-rawGal(1,1))/3600/24;
timGal=decimate(rawGal(:,1),floor(galSampF/sampF))';
inGal=(decimate(rawGal(:,2),floor(galSampF/sampF)));
outGal=decimate(rawGal(:,3),floor(galSampF/sampF));

% timGal=rawGal(:,1)';
% inGal=rawGal(:,2);
% outGal=rawGal(:,3);

geoGal = galR*[sin((90-galDec)*pi/180)*cos((galRA)*pi/180) sin((90-galDec)*pi/180)*sin((galRA)*pi/180) cos((90-galDec)*pi/180)];
trans = geoGal + [0 0 eaR];
labGal = trans*[cos((90-seaLat)*pi/180) 0 sin((90-seaLat)*pi/180);
    0 1 0;
    -sin((90-seaLat)*pi/180) 0 cos((90-seaLat)*pi/180)];


inGalAn = [];
for index = 1:length(timGal) 
    labGalAn = labGal*[cos(2*pi*timGal(index)-seaLong*pi/180) sin(2*pi*timGal(index)-seaLong*pi/180) 0;
    -sin(2*pi*timGal(index)-seaLong*pi/180) cos(2*pi*timGal(index)-seaLong*pi/180) 0;
    0 0 1];
    inGalAn = [inGalAn; -labGalAn(1)/galR];
end

figure(1)
l=plot(timGal,(inGal), timGal, (inGalAn));
ylabel('Basis function','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
legend('Astropy','Analytical')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
grid on
