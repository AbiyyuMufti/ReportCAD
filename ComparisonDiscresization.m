clear;
Tq = 1e-2;
SimTime = 0.1;
out = sim("AproxDiscreteComparison.slx");
subplot(1,3,1);
plot(out.tout, out.simout(1:end,1), out.tout, out.simout(1:end,2));
title("Forward Euler");
legend('Integrator', 'Forward Euler');
grid on;
subplot(1,3,2);
plot(out.tout, out.simout(1:end,1), out.tout, out.simout(1:end,3));
title("Backward Euler")
legend('Integrator', 'Backward Euler');
grid on;
subplot(1,3,3);
plot(out.tout, out.simout(1:end,1), out.tout, out.simout(1:end,4));
title("Trapezoidal (Tustin)")
legend('Integrator', 'Tustin');
grid on;