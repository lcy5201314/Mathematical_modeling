% ������һ������㷨������
% num:    ÿ���������Ŀ
% rClass: ����ֵ��x��Xr������ڵ��������
% xClass: ����ֵ��x���������
% ========================================================
function [rClass,xClass]=NNforCondense(Xr,x)

tic
% X = [randn(200,2)+ones(200,2);...
%      randn(200,2)-2*ones(200,2);...
%      randn(200,2)+4*ones(200,2);];
%  x=randn(1,2);%��������
[row,col]=size(Xr);
Xdist=zeros(row,1);
for i=1:row
    Xdist(i)=norm(x(1,1:2)-Xr(i,1:2))^2;
end
[Xdist,ind]=sort(Xdist,'ascend');
B=dist(1);
Xnn=Xr(ind(1),:);
rClass=Xnn(1,3);
xClass=x(1,3);
times=toc;