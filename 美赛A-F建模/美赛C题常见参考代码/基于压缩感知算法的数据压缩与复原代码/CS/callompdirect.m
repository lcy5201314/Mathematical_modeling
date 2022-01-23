clc;
clear;


[Y,fs,bits]=wavread('wrd.wav');

LL=85;
 M=30:40:280;

for n=1:7
   n
    [x1,x2]=ompdirect(M(n),Y);
    r(n)=x1;
    snr(n)=x2;

end
hold on;
figure(2)
plot(r,snr,'-o')
title('omp directly')