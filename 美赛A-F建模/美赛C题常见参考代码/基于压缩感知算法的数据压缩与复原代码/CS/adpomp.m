%novoiced and voiced decision omp ...
%this is a adaptive  cs-algorithm in speech signal...
%这个想法失败

function [r,snr]=adpomp(s1,s2)

N=320;
load('ZER.mat');
[Y,fs,bits]=WAVREAD('wrd.wav');
L=length(Y);
LL=length(Y)/N;
Yr=zeros(L,1);
Ys=abs(fft(Y));

%unvoiced and voiced decision
%then do omp

MM=zeros(1,LL);num=1;
lastpos=0;
    ii=1;
    jj=1 ;
    kk=1;
    Z1=zeros(N,1);
    Z2=zeros(N,1);
    B=zeros(L,1);
 while ii<=(LL-1)*N
 
       if Em(ii)<0.005                             %自己设定能量的阈值
          B(ii)=0; Yr(ii)=0;ii=ii+1;
     

       
       
       
       else
            if  ZER(ii) <=150                       %进行点数N1的 omp
%       for n=1:N
%             B(ii)=1; Z1(n)=Y(ii);ii=ii+1;
%       end
           B(ii:ii+N-1)=1;Z1=Y(ii:ii+N-1);ii=ii+N;
                  
          lastpos=ii;   Y1=Z1   ;               %  Z1在此处是行向量
                  
K=1; 
Ys=abs(fft(Y1));
for n=1:N
    
if  Ys(n)> s1                                     %设置浊音稀疏度门限
    K=K+1;
end
end

M=K*round(log2(N/K)+1)+30;
MM(num)=M;
num=num+1;
m=M/4;
% 恢复信号和原始信号对比
% m=12;
% M=80;
hat_x=omp(Y1,M,N,m);

          
          pos=lastpos-N;
                  for n=1:N
                    Yr(pos)=hat_x(n);
                    pos=pos+1;
                    
                  end 
           
             
            
           
        
            else                         %进行点数N2 omp
          
%          for n=1:N
%         B(ii)=-1;Z2(n)=Y(ii);ii=ii+1;
%          end
%            
         
          B(ii:ii+N-1)=-1;Z2=Y(ii:ii+N-1);ii=ii+N;
       lastpos=ii; Y2=Z2;pos=lastpos-N;
                
                 
             
K=1;
Ys=abs(fft(Y2));

for n=1:N
    
if  Ys(n)>s2                             %设置清音稀疏度门限
    K=K+1;
end
end

M=K*round(log2(N/K)+1)+10;
MM(num)=M;
num=num+1;
m=M/4;
% 恢复信号和原始信号对比
% m=12;
% M=40;
hat_x=omp(Y2,M,N,m);
 
             
             
             
              
                   for n=1:N
               Yr(pos)=hat_x(n);
               pos=pos+1;
                   end  
           
                 end
        
           
            
          
        
            end
    
 end
 
%  figure(2)
%  subplot(414);
%  plot(B,'.');%绘制浊音清音的鉴别区间
%  
%  
%  
%  
%  figure(1)
%  subplot(313)
%  plot(Yr);
%  
%  legend('The recovery speech');
% figure(3)
% plot(Y,'r');
% hold on;
% plot(Yr,'k.-')
% wavplay(Yr)
error=norm(Y-Yr)/norm(Y)
snr=10*log10(1/error);
r=sum(MM)/L;

         
        
