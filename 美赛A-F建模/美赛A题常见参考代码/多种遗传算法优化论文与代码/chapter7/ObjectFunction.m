function obj=ObjectFunction(X)
%% ���Ż���Ŀ�꺯��
% X��ÿ��Ϊһ������
col=size(X,1);
for i=1:col
    obj(i,1)=21.5+X(i,1)*sin(4*pi*X(i,1))+X(i,2)*sin(20*pi*X(i,2));
end
