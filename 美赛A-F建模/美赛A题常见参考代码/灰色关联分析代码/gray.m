clc;
close;
clear all;
x=xlsread('C:\\Users\\lenovo\\Desktop\\gray_data1.xlsx');
x=x(:,2:end)';
column_num=size(x,2);
index_num=size(x,1);

% 1�����ݾ�ֵ������
x_mean=mean(x,2);
for i = 1:index_num
    x(i,:) = x(i,:)/x_mean(i,1);
end
% 2����ȡ�ο����кͱȽ϶���
ck=x(1,:)
cp=x(2:end,:)
cp_index_num=size(cp,1);

%�Ƚ϶�����ο��������
for j = 1:cp_index_num
    t(j,:)=cp(j,:)-ck;
end
%���������С��
mmax=max(max(abs(t)))
mmin=min(min(abs(t)))
rho=0.5;
%3�������ϵ��
ksi=((mmin+rho*mmax)./(abs(t)+rho*mmax))

%4���������
ksi_column_num=size(ksi,2);
r=sum(ksi,2)/ksi_column_num;

%5�����������򣬵õ����r3>r2>r1
[rs,rind]=sort(r,'descend');
rs,rind
