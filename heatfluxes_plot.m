clear all
close all
clc

T = [386.06 332.70 321.16 316.29];
q_aero = [549.60 1.0992e3 1.0992e3 1.0992e3];
q_kin = [31.25 31.25 31.25 31.25];
q_drop = [ 6.1774e+03  4.9599e+03 4.4460e+03  4.1253e+03 ];
q_evap = [1.4537e+05 8.1039e+04 6.2026e+04 5.2095e+04 ];
q_conv = [2.9277e+04  4.7013e+04  4.2142e+04 3.9102e+04];

plot(T,q_aero, T,q_kin, T, q_drop, T,q_evap, T,q_conv, 'linewidth', 1.5)
set(gca, 'XDir','reverse')
legend('q_{aero}', 'q_{kin}','q_{drop}','q_{evap}','q_{conv}')
xlabel('T [K]'); ylabel('Heat Flux [W/m^2]')
grid on
