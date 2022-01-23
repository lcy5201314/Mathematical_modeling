
clc;
clear;
t=0:0.1:25;
b=98.7729*0.001;
fc=669.1851*0.001;
a=0.1;
gt=a*t.^3.*exp(-2*pi*b*t).*cos(2*pi*fc*t);
subplot(211); plot(t,gt);xlabel('t/ms');ylabel('幅度')
subplot(212);plot(t,abs(fft(gt)));xlabel('f/kHz');ylabel('幅度')%我只画了一个滤波器
