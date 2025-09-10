load("Uraw.mat")

[n,xh] = hist(Uraw,100);

x0 = [sum(n),7,2];
chi = @(x,xh)x(1)*((xh.*x(2)).^(x(3)/2-1).*exp(-(xh.*x(2))/2)/(2^(x(3)/2)*gamma(x(3)/2)));

[x,resnorm,~,exitflag,output] = lsqcurvefit(chi,x0,xh,n);

cdf = cumsum(chi(x,xh))/max(cumsum(chi(x,xh)));
cutoff = xh(min(find(cdf>0.99)))

figure(5)
bar(xh,n,'FaceAlpha',0.5)
hold on
plot(xh,chi(x,xh),'LineWidth',1.5)
plot([cutoff cutoff], [0 500],'k--','LineWidth',1.5)
hold off
xlabel('$\chi^2$','Interpreter', 'latex')
ylabel('Number','Interpreter', 'latex')
legend('Without Data Quality Cut','Fit','Cutoff','Interpreter', 'latex')
set(gca,'FontSize',16);
grid on

fig2=figure(5);
set(fig2,'Units','Inches');
pos = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig2,'EP_chiCut.pdf','-dpdf','-r1200')
