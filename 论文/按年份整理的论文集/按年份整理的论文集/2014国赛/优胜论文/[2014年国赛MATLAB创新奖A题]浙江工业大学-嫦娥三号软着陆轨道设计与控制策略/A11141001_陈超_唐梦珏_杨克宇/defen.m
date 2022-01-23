clc,clear;
load('q1.mat'); % ������ǰ��õĴֱ��ϵ�������Ϣ
load('w.mat');  % ������ǰ��õľ����ϵİ�ȫ�÷���Ϣ
%�������ϼ���÷�
mymax = 255;
mymin = 1;
for i=1:20
    for j=1:20
        if w(i,j)~=0
            w(i,j)=(w(i,j)-mymin)/(mymax-mymin);
        end
    end
end

d=zeros(20,20);
for i=1:20
    for j=1:20
        d(i,j)=sqrt((i-10)^2+(j-10)^2);
    end
end

dmax=max(max(d));
d=d/dmax;
for i=1:20
    for j=1:20
        if w(i,j)~=0
            f(i,j)=0.7*w(i,j)+0.3*d(i,j);
        else
            f(i,j)=1111;
        end
    end
end
x=30:50:980;
y=30:50:980;
rank=zeros(10,3);
for i=1:10
    [rank(i,1),rank(i,2),rank(i,3)]=myMin(f);
    f(rank(i,1),rank(i,2))=2222;
    rank(i,1)=x(rank(i,1))/10;% ��λת��
    rank(i,2)=y(rank(i,2))/10;% ��λת��
end

% �ֱ��ϼ���÷�
% mymax = 11.3696;
% mymin = 0.0247;
% for i=1:230
%     for j=1:230
%         if q1(i,j)~=-1
%             q1(i,j)=(q1(i,j)-mymin)/(mymax-mymin);
%         end
%     end
% end
% 
% d=zeros(230,230);
% for i=1:230
%     for j=1:230
%         d(i,j)=sqrt((i-115)^2+(j-115)^2);
%     end
% end
% 
% dmax=max(max(d));
% d=d/dmax;
% for i=1:230
%     for j=1:230
%         if q1(i,j)~=-1
%             f(i,j)=0.7*q1(i,j)+0.3*d(i,j);
%         else
%             f(i,j)=1111;
%         end
%     end
% end
% x=5:10:2295;
% y=5:10:2295;
% rank=zeros(10,3);
% for i=1:10
%     [rank(i,1),rank(i,2),rank(i,3)]=myMin(f);
%     f(rank(i,1),rank(i,2))=2222;
%     rank(i,1)=x(rank(i,1));
%     rank(i,2)=y(rank(i,2));
% end

