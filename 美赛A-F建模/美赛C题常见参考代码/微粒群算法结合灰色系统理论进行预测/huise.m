function [GM f]=huise(data,p)
N=1;%Ԥ�������
T=length(data);
data=data*p(2);
X0(1)=(3*data(1)+data(2))/4; %����ƽ������
for i=2:T-1
   X0(i)=(data(i-1)+2*data(i)+data(i+1))/4;
end
X0(T)=((data(T-1)+3*data(T))/4);
X1=(cumsum(X0,2))';%��ƽ��������һ���ۼ�
for i=1:T-1
    B(1,i)=-(X0(i+1)+X0(i))/2;
    B(2,i)=-((1-p(1))*X1(i+1)+p(1)*X1(i));   %����p(1)
    B(3,i)=1;
end
B=B';
for i=1:T-1
    Y(i,1)=(X0(i+1)-X0(i))/2;
end
H=B\Y;%����ϵ��
if (H(1)^2-4*H(2))>=0
 v1=(-H(1)+sqrt(H(1)^2-4*H(2)))/2;
 v2=(-H(1)-sqrt(H(1)^2-4*H(2)))/2;
 c2=(X0(2)+v1*exp(v1)*H(3)/H(2)-v1*X0(1))/(v2-v1);
 c1=X0(1)-(H(3)/H(2))-c2; 
else
  W=sqrt(4*H(2)-H(1)^2)/2;
 v1=-H(1)/2+W;
 v2=-H(1)/2-W;
 a1=X0(1)-H(3)/H(2);
 a2=(X0(2)+a1*H(1)/2)/W; 
end
    %��÷���ϵ��
if (H(1)^2-4*H(2))>=0
    for i=1:T+N
       G(i+1)=c1*exp(v1*i)+c2*exp(v2*i)+H(3)/H(2);
%       F(i)=X1(i+1)-G(i+1);
    end
else
    for i=1:T+N
       G(i+1)=exp(-H(1)*i/2)*(a1*cos(W*i)+a2*sin(W*i))+H(3)/H(2);
%        F(i)=X1(i+1)-G(i+1);
    end
end   

yu(1)=G(2)-data(1);
q(1)=data(1)-yu(1);
for i=2:T+N
       yu(i)=(G(i+1)-G(i))/p(2);  %Ԥ��ֵ
%        q(i)=data(i)-yu(i); %Ԥ��ֵ��ԭֵ�Ĳв�
%        p(i)=q(i)/data(i+1)*100; %������
end
data=data/p(2);
GM=minf(yu,data);%minf����ֵΪԤ��ֵK��ԭֵdata֮�������ۼ�ֵ
f=yu;





