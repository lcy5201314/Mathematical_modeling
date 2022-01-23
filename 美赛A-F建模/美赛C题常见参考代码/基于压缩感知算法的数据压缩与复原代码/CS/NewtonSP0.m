%牛顿梯度法
% function x_recovery=newton(M,N,x,m)
clc;
clear

K=6;      %  稀疏度(做FFT可以看出来)
N=320;    %  信号长度
M=100;     %  测量数(M>=K*log(N/K),至少40,但有出错的概率)

% m=10;     %迭代次数
% f1=50;    %  信号频率1
% f2=100;   %  信号频率2
% f3=200;   %  信号频率3
% f4=400;   %  信号频率4
% fs=800;   %  采样频率
% ts=1/fs;  %  采样间隔
% Ts=1:N;   %  采样序列
% x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  完整信号
% x=x';
[x,fs,bits]=WAVREAD('M1_1',[8400 8719]);
A=randn(M,N);
y=A*x;
for kk=2:N
    for nn=1:N
        dctbasis(kk,nn)=(2/N)^0.5*cos((2*(nn-1)+1)*(kk-1)*pi/2/N);
    end
end
for nn=1:N
    dctbasis(1,nn)=(1/N)^0.5*cos((2*(nn-1)+1)*(1-1)*pi/2/N);
end

B=dctbasis;                       
T=A*B'; %  傅里叶正变换矩阵 %y=A*x=A*B'*s=T*s;    s=B*x;
Aug_t=[];
index=[];
% Id=[];
rn=y;
% a=1;%步长

hat_x=zeros(N,1);
hat_xk=zeros(K,1);
hat_xs=zeros(N,1);
% G=Aug_t'*Aug_t;
% b=Aug_t'*y;
% f=1/2*hat_x'*G*hat_x-b'*hat_x;
m=5;
% while(norm(rn)>10^-2)
    for times=1:m
    for col=1:N
        inner(col)=abs(T(:,col)'*rn);
    end
    [val,pos] = sort(inner, 'descend');
 
    index=union(index,pos(1:K));
 
    
    Aug_t=T(:,index) ; %M*i
    
    l=size(Aug_t,2);
       
    

     H=-Aug_t'*Aug_t;                              %海森矩阵i*i
     gn=Aug_t'*rn;                                 %i*1,梯度
     d=-H^(-1)*gn;                                      %步长i*1，d=(rn'*cn)/norm(cn).^2;

     hat_xs(index)=hat_xs(index)+d  ;         %hat_x(pos_array)=hat_x(pos_array)+d*gn  
     hat_x=hat_xs;

   hat_x(hat_x==0)=[];
  [val,pos] = sort(abs(hat_x), 'descend');

  index_new=index(pos(1:K));


           
  
   rn=rn-Aug_t*d;                                  %a*Aug_t*d
%      T(:,index)=zeros(M,size(index,2)); 
% if norm(rn)<10^-3
%      break;
%  end
times=times+1;

    end
% Aug_t=T(:,index_new);
% xs=pinv(Aug_t)*y;
% hat_xs(index_new) =xs;
x_recovery=real(B'*hat_xs);
times
rn=norm(rn)
error=norm(x_recovery-x)^2/norm(x)^2  ;
snr=10*log10(1/error)
plot(x_recovery,'k.-')
hold on;
plot(x,'r');
legend('x_recovery','x');
hold off