clear; clc;

% Creation of model's parameters
x = [0.2244 0.5636 16.6338];
k1 = x(1);
k2 = x(2);
C12 = x(3);

syms s z Kd Kp Ki Tq

s = ((z-1)*2)/((z+1)*Tq);
HCs = collect((Kd*s^2 + Kp*s + Ki)/(Tq*s^2 + s), 'z');
% pretty(HCs)
Pre_fill = collect(Ki/(Kd*s^2 + Kp*s + Ki), 'z');
% pretty(Pre_fill);
% [num_, denum_] = numden(Pre_fill);
% A_cc = eval(coeffs(denum_,z))
% B_cc = eval(coeffs(num_,z))

filter = false;
% Coefficients for continous time PID controler
if filter 
    object = 'DigitalPIDWithFilter.slx';
    Kp = 33.4244;
    Ki = 114.9607;
    Kd = 3.1197;
    
else
    object = 'DigitalPID.slx';
    Kp = 3;
    Ki = 10;
    Kd = 1;
end
Tq   = 0.001;

nz2 = (Ki*Tq^2 + 2*Kp*Tq + 4*Kd);
nz1 = (2*Ki*Tq^2 - 8*Kd);
nz0 = Ki*Tq^2 - 2*Kp*Tq + 4*Kd;

dz2 = 6*Tq;
dz1 = -8*Tq;
dz0 = 2*Tq;

n_prz0 = Ki*Tq^2;
n_prz1 = 2*Ki*Tq^2;
n_prz2 = Ki*Tq^2;

d_prz0 = Ki*Tq^2 - 2*Kp*Tq + 4*Kd;
d_prz1 = 2*Ki*Tq^2 - 8*Kd;
d_prz2 = Ki*Tq^2 + 2*Kp*Tq + 4*Kd;


% start simulation of the DC motor's model
SimTime = 1;
Ts = 1e-3;
out = sim(object);
plot(out.tout, out.simout(1:end,1)); hold on; 
plot(out.tout, out.simout(1:end,2), 'r'); hold on;
plot(out.tout, out.simout(1:end,3), 'g--'); grid on;
legend('Continuous', 'DiscretFilter', 'With Function');
