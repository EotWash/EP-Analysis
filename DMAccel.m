%% Acceleration towards DM
G = 6.67430e-11; % Gravitational constant (m^3/kg/s^2)
Msun = 1.9891e30; % Mass of sun (kg)
kpc = 3.09e19;

%% Einasto Profile from https://arxiv.org/abs/2303.12838

M = 0.62e11*Msun;
beta = 0.91;
R = 3.86*kpc;

x1 = 7.5*kpc;
x2 = 8.5*kpc;
fun = @(r) r.*exp(-(r/R).^beta);
grad = (integral(fun, 0, x2)-integral(fun, 0, x1))/(x2-x1);
(G*M/R^3*grad)

%% Einasto Profile from https://arxiv.org/pdf/1906.06133

M2 = 2.8e11*Msun;
alpha = 0.18;
R2 = 182*kpc;
Rs = 11*kpc;

fun = @(r) r.^2.*exp(-(2/alpha).*((r/Rs).^alpha-1));
rho0 = M2/(4*pi*integral(fun, 0, R2));

x1 = 7.5*kpc;
x2 = 8.5*kpc;
fun = @(r) r.*exp(-(2/alpha)*((r/Rs).^alpha-1));
grad = (integral(fun, 0, x2)-integral(fun, 0, x1))/(x2-x1);
4*pi*G*rho0*grad

%% gNFW Profile from https://arxiv.org/pdf/1906.06133

M2 = 5.5e11*Msun;
gamma = 1.2;
R2 = 174*kpc;
Rs = 9*kpc;

fun = @(r) r.^2.*(Rs./r).^gamma.*(1+r/Rs).^(gamma-3);
rho0 = M2/(4*pi*integral(fun, 0, R2));

x1 = 7.5*kpc;
x2 = 8.5*kpc;
fun = @(r) r.*(Rs./r).^gamma.*(1+r/Rs).^(gamma-3);
grad = (integral(fun, 0, x2)-integral(fun, 0, x1))/(x2-x1);
4*pi*G*rho0*grad