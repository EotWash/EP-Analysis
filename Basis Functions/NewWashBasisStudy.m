sampF = 2.2866e-04;
fDay = 1/24/3600;
wDay = 2*pi*fDay;

if true
    % Loading in galaxy basis funtions outputted from galVect.py
    rawGal=load('galVectMinTest.out');
end
galSampF = 1/(rawGal(2,1)-rawGal(1,1))/3600/24;
timGal=rawGal(:,1);
inGal=rawGal(:,2);
outGal=rawGal(:,3);

% timGal = timGal - (timGal>307.042)/24 - (timGal<69.082)/24 - (timGal<433.041)/24 - (timGal>671.041)/24;

timGalD=decimate(timGal,floor(galSampF/sampF));
inGalD=(decimate(inGal,floor(galSampF/sampF)));
outGalD=(decimate(outGal,floor(galSampF/sampF)));


figure(1)
l=plot(timGal,inGal);%, timGalD, inGalD);
ylabel('Basis function','Interpreter', 'latex')
xlabel('Time (days)','Interpreter', 'latex')
legend('Raw','Decimated')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
grid on

figure(2)
l=plot(timGal,inGal-max(inGal)*sin(2*pi*timGal), timGalD, inGalD-max(inGal)*sin(2*pi*timGalD));
ylabel('Residual','Interpreter', 'latex')
xlabel('Time (s)','Interpreter', 'latex')
legend('Raw','Decimated')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
grid on