etaSun = 2.8e-13;
delta_BL = 0.0359;
he_BL = 1-2/4;
% sun_BL = 0.28*he_BL;
sun_BL = 1;
Rsun = 149.6e9; % Radius from earth to sun (m)

etaPast = 1.3e-13*2;
REarth = 6378e3;

etaMicro = 2.7459e-15*2;
RMicro = 712e3;
fe_BL = 1-26/55;
o_BL = 1-8/16;
si_BL = 1-14/28;
mg_BL = 1-12/24.3;
s_BL = 1-16/32;
ni_BL = 1-28/58.7;
ca_BL = 1-20/40;
al_BL = 1-13/27;
% earth_BL = 0.321*fe_BL+0.301*o_BL+0.151*si_BL+0.139*mg_BL...
%     +0.029*s_BL+0.018*ni_BL+0.015*ca_BL+0.014*al_BL;
earth_BL = 1;
micro_BL = 0.0583;

lambda = logspace(9.75, 12);

alpha = etaSun./(delta_BL*sun_BL*(1+Rsun./lambda).*exp(-Rsun./lambda));
alphaM = etaMicro./(micro_BL*earth_BL*(1+RMicro./lambda).*exp(-RMicro./lambda));
alphaP = etaPast./(delta_BL*earth_BL*(1+REarth./lambda).*exp(-REarth./lambda));


figure(11)
l=loglog(lambda, alpha, lambda, alphaP, lambda, alphaM);
xlabel('$\lambda$ (m)','Interpreter', 'latex')
ylabel('$\alpha$','Interpreter', 'latex')
set(gca,'FontSize',16);
set(l,'LineWidth',1.5);
grid on
ylim([1e-15 1e-5])
legend('Current Results', 'Past Earth Results', 'MICROSCOPE','Interpreter', 'latex')

