clear; clc;
disp('Busy!');

number = 100;
Tq = 0.1;
t = 0:Tq:10-Tq;
points_number = size(t,2);
P=zeros(points_number, number);

xInit = [0.5 0.5 10];
fieldSize = 0.7;

min = xInit*(1-fieldSize);
max = xInit*(1+fieldSize);

T = min.*ones(number, 3) + (max - min).*rand(number, 3);
SimTime = 10;
Ts = 1e-3;
for i=1:number
    k1 = T(i, 1);
    k2 = T(i, 2);
    C12 = T(i, 3);
    [num, denum] = linmod('DoubleMassSimplified');
    sys = tf(num, denum);
    P(:,i) = step(sys,t);
end

plot(P);
grid on;
save reference P T points_number Tq

clc;
disp('Done!');
