clc;
clear;
close all;
M=70;
N=320;
m=62;
[x,fs,bits]=WAVREAD('M1_1',[8400 8719]);
plot(x,'r')
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            newton
% x_recovery0=newton(M,N,x,m);                           
% 
% errornewton=norm(x-x_recovery0).^2/norm(x).^2
% snr=10*log10(1/errornewton)
% plot(x_recovery0,'k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                steepest
x_recovery1=steepest(M,N,x,m);
errorsteepest=norm(x-x_recovery1).^2/norm(x).^2

snrsteepest=10*log10(1/errorsteepest)
plot(x_recovery1,'g')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 omp
% x_recovery2=omp(M,N,x,m);
% erroromp=norm(x-x_recovery2).^2/norm(x).^2
% 
% snromp=10*log10(1/erroromp)
% 
% plot(x_recovery2,'b')











legend('orignal ','newton','steepest','omp');
hold off;
