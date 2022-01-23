
% function [r,snr]=decisionomp(Y,M1,M2)
clc;clear
[Y,fs,bits]=wavread('wrd.wav');
N=320;
% [Y,fs,bits]=wavread('wrd.wav');
% plot(Y)
L=length(Y);



Fr=enframe(Y,N);
LL=size(Fr,1);
MM=zeros(1,LL);
x=zeros(1,N);
hat_x=zeros(LL,N);
k=6;
M1=80;
M2=30;

for ii=1:LL
    x=Fr(ii,:);
    
    
    
    
    %compute Zn,Em
    
   Em=sum(x.*x);
   s=0;
   for n=1:N-1
   s=s+abs(sign(x(n+1)-sign(x(n))));
   end 
   Zn=0.5*s;
   ZER=Zn/Em;
   
   
if Em<0.005
    hat_x(ii,:)=zeros(1,N);
    MM(ii)=0;
    
else
    
    
    if ZER<150
        m=M1/4;
  hat_x(ii,:)=omp(x',M1,N,m);
  MM(ii)=M1;
  
   
    else
        
        
        m=M2/4;
        
      hat_x(ii,:)=omp(x',M2,N,m); 
      MM(ii)=M2;
    end
        
end



end

        
        
        
        
        
        
        
        
%   x1=abs(fft(x));
%     k=1;
%     for jj=1:N
%         if x1>0.2
%             k=k+1;
%         end
%     end

%     M=k*round(log2(N/k))+20;
  
    

Y_recovery=overlapadd(hat_x);


% 
% plot(Y_recovery,'r.-');
% legend('orignal',' recovery')

%plot(abs(fft(Fr(20,:))));
% wavplay(Y_recovery);


error=norm(Y-Y_recovery)/norm(Y)

r=(sum(MM)/LL)/N
snr=10*log10(1/error);
wavwrite(Y_recovery,16000,'2');
    