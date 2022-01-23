T1=0.5:0.5:4;
T2=0.05:0.05:0.4;
clc;


for n=1:7
    n
    [x1,x2]=adpomp(T1(n),T2(n));
    r(n)=x1;
    snr(n)=x2;
end
  close all;  
plot(r,snr,'o-');
