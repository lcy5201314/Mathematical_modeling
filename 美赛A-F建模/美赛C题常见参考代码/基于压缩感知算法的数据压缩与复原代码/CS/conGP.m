%共轭梯度算法


K=7;      %  稀疏度(做FFT可以看出来)
N=256;    %  信号长度
M=64;     %  测量数(M>=K*log(N/K),至少40,但有出错的概率)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%m

for times=1:m
    for col=1:N
        inner(col)=abs(T(:,col)'*rn);
    end
    [val,pos]=max(inner);
    Aug_t=[Aug_t,T(:,pos)];  %M*i
    pos_array(times)=pos; 
    T(:,pos)=zeros(M,1); 
    
    
    
    
  
    gn=Aug_t'*rn;   %i*1,梯度
    
    cn=Aug_t*gn;    %M*1
    ak=(rn'*cn)/norm(cn).^2;%步长
    
    
    hat_x(pos_array)=hat_x(pos_array)+d*gn  ;
   
    rn=rn-d*Aug_t*gn;
end

x_recovery=real(B'*hat_x);
    
    






