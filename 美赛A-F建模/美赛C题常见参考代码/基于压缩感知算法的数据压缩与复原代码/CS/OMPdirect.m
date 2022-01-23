% close all
% clc;
% % clear;

% function [r,snr]=ompdirect(M,Y)

[Y,fs,bits]=wavread('wrd.wav');
N=320;
M=128;

L=length(Y);

Fr=enframe(Y,N);
LL=size(Fr,1);



x=zeros(1,N);
hat_x=zeros(LL,N);


m=M/4;
for ii=1:LL
    x=Fr(ii,:);
    
    
%     
%     
%     compute Zn,Em
%     
%    Em=sum(x.*x);
%    s=0;
%    for m=1:N-1
%    s=s+abs(sign(x(m+1)-sign(x(m))));
%    end 
%    Zn=0.5*s;
%    ZER=Zn/Em;
%    
%    
% if Em<0.036
%     hat_x(ii,:)=zeros(1,N);
%     
% else
%     if ZER<100
%   hat_x(ii,:)=omp(x',M,N,m);
%     else
%       hat_x(ii,:)=omp(x',M,N,m); 
%     end
%         
% end


hat_x(ii,:)=omp(x',M,N,m);
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

% plot(Y,'r')
% hold on;
% plot(Y_recovery,'k.-');
% legend('orignal',' recovery')

%plot(abs(fft(Fr(20,:))));



error=norm(Y-Y_recovery)/norm(Y)

r=M/N;
snr=10*log10(1/error);
wavwrite(Y_recovery,16000,'1')

    