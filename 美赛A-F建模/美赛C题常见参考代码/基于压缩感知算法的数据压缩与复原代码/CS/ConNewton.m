%�����ݶ�


K=7;      %  ϡ���(��FFT���Կ�����)
N=256;    %  �źų���
M=64;     %  ������(M>=K*log(N/K),����40,���г���ĸ���)
f1=50;    %  �ź�Ƶ��1
f2=100;   %  �ź�Ƶ��2
f3=200;   %  �ź�Ƶ��3
f4=400;   %  �ź�Ƶ��4
fs=800;   %  ����Ƶ��
ts=1/fs;  %  �������
Ts=1:N;   %  ��������
x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts);  %  �����ź�
x=x';
A=randn(M,N);
y=A*x;
B=fft(eye(N,N))/sqrt(N); %y=A*x=A*B'*s=T*s;    s=B*x;
T=A*B';
Aug_t=[];
Id=[];
rn=y;


hat_x=zeros(N,1);


    for col=1:N
        inner(col)=abs(T(:,col)'*rn);
    end
    [val,pos]=max(inner);
    
    Aug_t=[Aug_t,T(:,pos)];  %M*i
    G=Aug_t'*Aug_t;
    pos_array(times)=pos; 
    T(:,pos)=zeros(M,1);  
    Id=[Id,pos];
    g1=Aug_t'*rn;
    
    
 
    s1=g1;%i*1,���ݶ�
    a1=g1'*s1/2*s1'*G*s1;   %����
    hat_x(pos_array)=hat_x(pos_array)+a1*s1  ;
 for iter=2:20
         for col=1:N
        inner(col)=abs(T(:,col)'*rn);
         end
    [val,pos]=max(inner);
    
    Aug_t=[Aug_t,T(:,pos)];  %M*i
    new
    G=Aug_t'*Aug_t;
    pos_array(times)=pos; 
    T(:,pos)=zeros(M,1);  
         gn=Aug_t'*rn;
         sk=gn+norm()
         ak=gk'*sk/2*sk'*G*sk;
             hat_x(pos_array)=hat_x(pos_array)+ak*sk  ;
    
   
    rn=rn-as*Aug_t*sk;




 end

x_recovery=real(B'*hat_x);
plot(x_recovery,'k.-');
hold on;
plot(x,'r');



legend('x','x_recovery')
    
    
        toc




