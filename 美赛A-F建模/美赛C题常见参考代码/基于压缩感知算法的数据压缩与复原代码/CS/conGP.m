%�����ݶ��㷨


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%m

for times=1:m
    for col=1:N
        inner(col)=abs(T(:,col)'*rn);
    end
    [val,pos]=max(inner);
    Aug_t=[Aug_t,T(:,pos)];  %M*i
    pos_array(times)=pos; 
    T(:,pos)=zeros(M,1); 
    
    
    
    
  
    gn=Aug_t'*rn;   %i*1,�ݶ�
    
    cn=Aug_t*gn;    %M*1
    ak=(rn'*cn)/norm(cn).^2;%����
    
    
    hat_x(pos_array)=hat_x(pos_array)+d*gn  ;
   
    rn=rn-d*Aug_t*gn;
end

x_recovery=real(B'*hat_x);
    
    






