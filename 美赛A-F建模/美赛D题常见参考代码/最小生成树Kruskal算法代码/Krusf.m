function [T c]=Krusf(d,flag)
% function [T c]=Krusf(d,flag)
% ����С��������Kruskal�㷨
%  ��Ȩ����Ĳ���������
%  1)һ��ı�Ȩ����Ϊnxnά�����÷�ʽ[T c]=Krusf(d)
%  2����Ȩ�����ǰ���зֱ��¼ͼ�����бߵ���ʼ�������ֹ���㣬
%     ����߲��ظ���¼�������м�¼��Ӧ�ߵ�Ȩֵ�����÷�ʽΪ[T c]=Krusf(d,1)
%  c:�������ķ���;
%  T:�������ı߼���;

if nargin==1
n=size(d,2);
m=sum(sum(d~=0))/2;
b=zeros(3,m);
k=1;
for i=1:n
    for j=(i+1):n
        if d(i,j)~=0
            b(1,k)=i;b(2,k)=j;
            b(3,k)=d(i,j);k=k+1;
        end
    end
end
else
    b=d;
end

n=max(max(b(1:2,:)));
m=size(b,2);
[B,i]=sortrows(b',3);
B=B';
c=0;T=[];
k=1;t=1:n;
for i=1:m
    if t(B(1,i))~=t(B(2,i))
         T(1:2,k)=B(1:2,i);
         c=c+B(3,i);
         k=k+1;
         tmin=min(t(B(1,i)),t(B(2,i)));
         tmax=max(t(B(1,i)),t(B(2,i)));
         for j=1:n
             if t(j)==tmax
                 t(j)=tmin;
             end
         end
    end
    if k==n
        break;
    end
end
T;
c;