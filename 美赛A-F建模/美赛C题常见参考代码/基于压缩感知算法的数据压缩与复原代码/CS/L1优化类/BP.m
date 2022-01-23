clear;clc;
m=400;k=10;
v1=zeros(1,m);
v=zeros(1,m);%��������ȫ��0��������
n=4*k;
t=0;
while  t<k
    u1=ceil((m-1)*rand(1,1));
    v1(u1)=1;
    t=sum(v1);
end
s=1;
for i=1:m
    if v1(i)==1
        S(s)=i; %x0��ϡ��ϵ����λ��
        s=s+1;
    end
end

x0=zeros(m,1);
for i=1:m
  for j=1:k
      if i==S(j)
         x0(i)=rand(1);%����һ�����ϡ���ź�
      end
  end
end


B=rand(n,m);%�����۲����B
y=B*x0;%measurement

%A*x=b x>=0 (LP ת��)
c=ones(2*m,1);

A=[B,-B];
b=y;
lb=ones(2*m,1);
options = optimset('LargeScale','off','Simplex','on');
[z,fval,exitflag,output] = linprog(c,[],[],A,b,lb,[],[],options);
output
u=z(1:m,1);
v=z(m+1:2*m,1);
x1=u-v;
E1=norm(x1-x0)

plot(x1,'k.-');hold on;
plot(x0,'r')
legend('recovery','orignal')















