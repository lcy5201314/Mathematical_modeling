% %内点惩罚牛顿法
% %min x1.^2+x2.^2    s.t 2*x1+x2-2<=0,x1>=1.
%X'=X-a*H'*gradient f(x),此程序a=1；
clc
m=zeros(1,50);a=zeros(1,50);b=zeros(1,50);f0=zeros(1,50);
syms x1 x2 e;
m(1)=1;c=0.2;a(1)=2;b(1)=-3;
f=x1^2+x2^2-e*(1/(2*x1+x2-2)+1/(1-x1)); f0(1)=15;
fx1=diff(f,'x1');fx2=diff(f,'x2');fx1x1=diff(fx1,'x1');fx1x2=diff(fx1,'x2');fx2x1=diff(fx2,'x1');fx2x2=diff(fx2,'x2');
for k=1:100
    x1=a(k);x2=b(k);e=m(k);
    for n=1:100
        
        f1=subs(fx1);
        f2=subs(fx2);
        f11=subs(fx1x1);
        f12=subs(fx1x2);
        f21=subs(fx2x1);
        f22=subs(fx2x2);
        if(double(sqrt(f1^2+f2^2))<=0.002)
            a(k+1)=double(x1);b(k+1)=double(x2);f0(k+1)=double(subs(f));
            break;
        else
            X=[x1 x2]'-inv([f11 f12;f21 f22])*[f1 f2]';
            x1=X(1,1);x2=X(2,1);
        end
    end
if(double(sqrt((a(k+1)-a(k))^2+(b(k+1)-b(k))^2))<=0.001)&&(double(abs((f0(k+1)-f0(k))/f0(k)))<=0.001)
      a(k+1)
      b(k+1)
      k
      f0(k+1)
      break;
else
      m(k+1)=c*m(k);
end
end