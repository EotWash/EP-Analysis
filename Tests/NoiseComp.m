len = length(longDat);

first = (longDat(1:len/2));
second = (longDat(len/2:len));

first = first(not(isnan(first)));
second = second(not(isnan(second)));

(std(second)/std(first)-1)*100


figure(99)
[NC,XC] = hist(first*1e15/(r*m),10);
[NS,XS] = hist(second*1e15/(r*m),XC);
bar(XC,(NC),'FaceAlpha',0.5)
hold on
bar(XS,(NS),'FaceAlpha',0.5)
hold off
% xlim([-4 4]*1e-25)
% xlabel('$g_{B-L}/\sqrt{\hbar c}$','Interpreter', 'latex')
ylabel('Number','Interpreter', 'latex')
% legend('Cosine','Sine','Thermal Noise','Interpreter', 'latex')
set(gca,'FontSize',16);
grid on
