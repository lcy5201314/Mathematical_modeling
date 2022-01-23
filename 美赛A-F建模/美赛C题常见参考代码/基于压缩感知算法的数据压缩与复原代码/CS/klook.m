clc;
clear;
close all;

[Y,fs,bits]=WAVREAD('wrd.wav');
N=320;


L=length(Y);


% Fr=enframe(Y,N);
% LL=size(Fr,1);
% x=zeros(1,N);
% 
% f=77;
    figure(1)
    x=Y(6400:6719);
   
    
     y2=abs(fft(x));
   
    plot(x)
      title('x')
    figure(2)

   
    plot(y2)
    axis([0 160 0 25]);
         title('fft') ;
         figure(3)
         plot(Y)