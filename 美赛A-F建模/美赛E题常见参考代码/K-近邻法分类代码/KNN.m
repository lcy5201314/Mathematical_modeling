% ====================K�����ڷ�(KNN)=================================
% X:  ѵ������
% x:  ��������
% K:  ������Ŀ 
% flag1: ��¼K������������ڵ�һ��ĸ���
% flag2: ��¼K������������ڵڶ���ĸ���
% ===================================================================
clear,close all;
N=150;
X = [randn(N,2)+2*ones(N,2);...
     randn(N,2)-2*ones(N,2);];
% ====================================================================
fig=figure;
plot(X(1:N,1),X(1:N,2),'r.')
hold on,plot(X(N+1:2*N,1),X(N+1:2*N,2),'b.')
title('��ʼ�����ֲ�ͼ'),axis([-10,10,-10,10])
% ====================================================================
x=randn(1,2);%��������
hold on,plot(x(1),x(2),'m+','MarkerSize',10,'LineWidth',2)
for i=1:2*N
    dist(i)=norm(x-X(i,:));
end
[Sdist,index]=sort(dist,'ascend');
K=5; %������Ŀ 
for i=1:K
    hold on,plot(X(index(i),1),X(index(i),2),'ko');
end
legend('Cluster 1','Cluster 2','x','Location','NW')
flag1=0;flag2=0;
for i=1:K
    if ceil(index(i)/N)==1
        flag1=flag1+1;
    elseif ceil(index(i)/N)==2
        flag2=flag2+1;
    end
end
disp(strcat('K�����а���',num2str(flag1),'����һ������'));
disp(strcat('K�����а���',num2str(flag2),'���ڶ�������'));
if flag1>flag2
    disp('�жϴ����������ڵ�һ��');
else
    disp('�жϴ����������ڵڶ���');
end
