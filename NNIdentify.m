clear; clc;
load("var4.mat");
load("nn3.mat");

P = Yout(1:100:end);
t = t(1:100:end);
k = 0.005; % noise
min_noise = -P(size(P,1))*k;
max_noise =  P(size(P,1))*k;
noise = min_noise + (max_noise - min_noise)*rand(size(P,1),1);
P_n = P+noise;
SimTime = 10;
Xnt=sim(net,P)

figure(1);
std_to_var = DoubleMassRunner(Xnt,true)
 
figure(2);
Xnt_noise=sim(net,P_n)
k1 = Xnt_noise(1);
k2 = Xnt_noise(2);
C12 = Xnt_noise(3);
val = sim('DoubleMassSimplified.slx');
val = val.simout;
val = val(1:100:end);
plot(t, P, 'b', t, P_n, 'r--', t, val, 'g.')
legend('Reference', 'With Noise', 'Identified');
% std_with_noise = DoubleMassRunner(Xnt_noise,true)

