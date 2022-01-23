clear;
clc
t=0:0.1:1;
x=3*sin(100*pi*t)+2*sin(68*pi*t)+5*sin(72*pi*t)+randn(1,length(t));

coefs=cwt(x,[1:0.2:3],'db3','plot');
title('对不同尺度小波系数')
ylabel('a')
xlabel('b')