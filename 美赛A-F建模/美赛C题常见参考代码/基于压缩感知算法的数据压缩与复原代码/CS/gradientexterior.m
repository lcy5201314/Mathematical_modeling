%����ݶȷ�(�����½��������ݶȷ����Ͻ���������Ѱ���½����Ĳ���)
%min f=x1.^2+x2.^2     s.t x1>=1;
clc
clear
m=zeros(1,50);a=zeros(1,50);b=zeros(1,50);f0=zeros(1,50);         
syms d x1 x2 e;                                                  
m(1)=1;c=10;a(1)=0;b(1)=0;                                        
f=x1^2+x2^2+e*(1-x1)^2;f0(1)=1;          %���ͷ����ʽ 
fx1=diff(f,'x1');                                                
fx2=diff(f,'x2');
for k=1:100                                                      
x1=a(k);x2=b(k);e=m(k);
for n=1:100                                                    
    f1=subs(fx1);
    f2=subs(fx2);%�����ݶ�
    if(double(sqrt(f1^2+f2^2))<=0.002)                          
        a(k+1)=double(x1);b(k+1)=double(x2);f0(k+1)=double(subs(f));
       break;
    else
       D=(x1-d*f1)^2+(x2-d*f2)^2+e*(1-(x1-d*f1))^2; %����һά������ʹ��Ѱ�Ҳ���d
       Dd=diff(D,'d'); dd=solve(Dd); x1=x1-dd*f1; x2=x2-dd*f2;%��һ���ŵ�
    end
end
if(double(sqrt((a(k+1)-a(k))^2+(b(k+1)-b(k))^2))<=0.001)&&(double(abs((f0(k+1)-f0(k))/f0(k)))<=0.001) %������ֹ��������һ�����ŵ������һ������룬������ֵ��������        
   a(k+1)
   b(k+1)
   k
   f0(k+1)
   break;
else
   m(k+1)=c*m(k);
end
end

