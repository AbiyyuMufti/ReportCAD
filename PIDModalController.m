clear;
close all;
clc;

k1 = 0.2244;
k2 = 0.5636;
C12 = 16.6338;
Ts = 1e-3;
SimTime = 5;
Omega = 10;

[num, denum] = linmod('DoubleMassSimplified');

% reduction
poles = sort(roots(denum))
polesQ = poles(1:end-1)
A_q = poly(polesQ);
B_q = A_q(end)*(num(end)/denum(end));

Ho = tf(num, denum)
Hv = tf(B_q, A_q)
t = [0:Ts:SimTime-Ts];
[orig.y, orig.t] = step(Ho,t);
[redu.y, redu.t] = step(Hv,t);

figure('Name','Reduction')
plot(orig.t, orig.y); grid on; hold on;
plot(redu.t, redu.y, 'r--'); grid on; hold off;
xlabel('Time (s)')
ylabel('Amplitude');
title('Transfer Function Reduction');
legend('Original','Reduction');

b0 = B_q(1);
a2 = A_q(1); 
a1 = A_q(2); 
a0 = A_q(3);

a3j = 1; % 
a2j = 3; % 
a1j = 3; % 
a0j = 1; % 



Kp = (a1j*Omega^2 - a0)/b0
Ki = (a0j*Omega^3)/b0
Kd = (a2j*Omega - a1)/b0

figure(2);
out = sim('PIDControlDMS');
plot(orig.t, orig.y, 'g'); grid on; hold on;
plot(out.tout, out.simout(1:end,1), 'r'); hold on;

out2 = sim('PIDControlDMSFiltered');
plot(out2.tout, out2.simout(1:end,1), 'b')
legend('Original','PID', 'PID Filtered');
xlabel('Time (s)')
ylabel('Amplitude');

