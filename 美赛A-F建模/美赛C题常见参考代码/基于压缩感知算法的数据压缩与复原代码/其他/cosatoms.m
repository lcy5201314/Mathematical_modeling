clear;
clc;

n=[1:256];
fx=[0:2*pi/10240:2*pi];

fm(1)=fx(45),fm(2)=fx(10),fm(3)=fx(24);

for i=1:3
f=cos(2*pi*fm(i)*n);

plot(f);
hold on;


end