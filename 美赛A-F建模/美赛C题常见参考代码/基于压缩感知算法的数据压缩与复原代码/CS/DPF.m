%最速度下降法
clc
clear
% function x_recovery=Steepest(M,N,x,m)
K=7;      %  稀疏度(做FFT可以看出来)
N=256;    %  信号长度
M=64;     %  测量数(M>=K*log(N/K),至少40,但有出错的概率)
m=2*K;
f1=50;    %  信号频率1
f2=100;   %  信号频率2
f3=200;   %  信号频率3
f4=400;   %  信号频率4
fs=800;   %  采样频率
ts=1/fs;  %  采样间隔
Ts=1:N;   %  采样序列
x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  完整信号
x=x';
A=randn(M,N);
y=A*x;
B=fft(eye(N,N))/sqrt(N); %y=A*x=A*B'*s=T*s;    s=B*x;
T=A*B';
Aug_t=[];
Id=[];
rn=y;


hat_x=zeros(N,1);
A=eye(N,N);
for times=1:m
for col=1:N
    inner(col)=abs(T(:,col)'*rn);
end
[val,pos]=max(inner);
Aug_t=[Aug_t,T(:,pos)];        %M*i
pos_array(times)=pos;
T(:,pos)=zeros(M,1);
if times==1
    g1=Aug_t'*rn;              %i*1,梯度
    go=g1;

    d1=-A*g1;
    a1=g1'*d1/d1'*A*d1;        %步长
    hat_x_old=hat_x;
    hat_x(pos_array)=hat_x(pos_array)+a1*d1 ;
    hat_x_new=hat_x;
                              %     deltax=hat_x_new-hat_x_old;
    rn=rn-a1*d1;
    A=eye(2);
else
    gn=Aug_t'*rn;
    if norm(gn)<10^3
        break;
    end
    deltag=gn-go;
    go=gn;
    
   
    dn=-A*gn;
    an=(gn'*dn)/(dn'*A*dn);
    hat_x(pos_array)=hat_x(pos_array)+an*dn ;
  
    
    Ek=((deltax*deltax')/(deltax'*deltag))-((A*deltag'*A)/(deltag'*A*deltag));
    A=A+Ek;
    rn=rn-an*dn;
end


end

x_recovery=real(B'*hat_x);
plot(x_recovery,'k.-');
hold on;
plot(x,'r');

legend('x','x_recovery')



