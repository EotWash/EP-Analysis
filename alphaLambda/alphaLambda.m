% Profile from https://arxiv.org/abs/2303.12838

G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
Msun = 1.9891e30; % Mass of sun (kg)
kpc = 3.09e19;

M = 0.62e11*Msun;
beta = 0.91;
R = 3.86*kpc;

a = 2*0.64e-15;

d = logspace(-3,3,100)*kpc;

l = 1e6*kpc;

fun = @(r) r.*exp(-r/l-(r/R).^beta);

I = [];

for in = d
    I = [I integral(fun, 0, in)];
end

dDiff = d(2:end);
iDiff = diff(I)./diff(d);

%%

lambda = logspace(-1,3,1000)*kpc;

x1 = 7.5*kpc;
x2 = 8.5*kpc;

alpha = [];

for in = lambda
    fun = @(r) r.*exp(-r/in-(r/R).^beta);
    grad = (integral(fun, 0, x2)-integral(fun, 0, x1))/(x2-x1);
    alphaI = a/(G*M/R^3*grad);
    alpha = [alpha alphaI];
end

%% Point Mass Approx
Rp = 8*kpc;
alphaP = a*Rp^2/(G*M)./(exp(-Rp./lambda).*(1+Rp./lambda));

%% Acceleration towards DM

x1 = 7.5*kpc;
x2 = 8.5*kpc;
fun = @(r) r.*exp(-(r/R).^beta);
grad = (integral(fun, 0, x2)-integral(fun, 0, x1))/(x2-x1);
(G*M/R^3*grad)


%%

figure(1)
l=loglog(d/kpc,I);
ylabel('Integral ($m^2$)','Interpreter', 'latex')
xlabel('Distance (kpc)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(l,'LineWidth',1.5);
grid on

figure(2)
l=loglog(dDiff/kpc,G*M/R^3*abs(iDiff));
ylabel('Acceleration ($m/s^2$)','Interpreter', 'latex')
xlabel('Distance (kpc)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(l,'LineWidth',1.5);
grid on

figure(4)
l=semilogx(d/kpc, M/4/pi/R^3*exp(-(d/R).^beta)/Msun*kpc^3);
ylabel('Density ($M_\odot/kpc^3$)','Interpreter', 'latex')
xlabel('Distance (kpc)','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(l,'LineWidth',1.5);
grid on


figure(3)
l=loglog(lambda/kpc, abs(alpha),lambda/kpc, abs(alphaP),[8 8],[1e-8 1e3],'k--');
hold on
text(7.2,0.6e-2,'Earth to Milky Way Distance','Interpreter', 'latex','FontSize',16,'Rotation',90)
hold off
ylabel('$\alpha$','Interpreter', 'latex')
xlabel('$\lambda$ (kpc)','Interpreter', 'latex')
legend('Dark Matter Halo','Point Mass Approx.','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'MarkerSize',16);
set(l,'LineWidth',1.5);
grid on
ylim([1e-7 1e2])
%%
if (false)
    fig2=figure(3);
    set(fig2,'Units','Inches');
    pos = get(fig2,'Position');
    set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(fig2,'EP_alphaLambda.pdf','-dpdf','-r1200')

end
