
%%ţ����㷨��min f=x1.^2+x2.^2     s.t x1>=1;


clc
m=zeros(1,50);a=zeros(1,50);b=zeros(1,50);f0=zeros(1,50);%a bΪ���ŵ����꣬f0Ϊ���ŵ㺯��ֵ��f1 f2���ŵ��ݶȡ�
syms x1 x2 e;                                            %eΪ�����ӡ�
m(1)=1;c=10;a(1)=0;b(1)=0;                               %cΪ����ϵ��������ֵ��
f=x1^2+x2^2+e*(x1-1)^2;f0(1)=1;
fx1=diff(f,'x1');fx2=diff(f,'x2');fx1x1=diff(fx1,'x1');fx1x2=diff(fx1,'x2');fx2x1=diff(fx2,'x1');fx2x2=diff(fx2,'x2');
for k=1:100                                                          %��㷨e����ѭ��.
    x1=a(k);x2=b(k);e=m(k);
    for n=1:100                                                      %�ݶȷ�������ֵ��
        f1=subs(fx1);                                                 %����ݶ�ֵ�ͺ�ɭ����
        f2=subs(fx2);
        f11=subs(fx1x1);
        f12=subs(fx1x2);
        f21=subs(fx2x1);
        f22=subs(fx2x2);
        if(double(sqrt(f1^2+f2^2))<=0.001)                                  %����ֵ��������
            a(k+1)=double(x1);b(k+1)=double(x2);f0(k+1)=double(subs(f));
            break;
        else
            X=[x1 x2]'-inv([f11 f12;f21 f22])*[f1 f2]';
            x1=X(1,1);x2=X(2,1);
        end
    end
if(double(sqrt((a(k+1)-a(k))^2+(b(k+1)-b(k))^2))<=0.001)&&(double(abs((f0(k+1)-f0(k))/f0(k)))<=0.001)   %�����ӵ�����������
      a(k+1)   %������ŵ����꣬�����ӵ�������������ֵ
      b(k+1)
      k
      f0(k+1)
      break;
else
      m(k+1)=c*m(k);
end
end