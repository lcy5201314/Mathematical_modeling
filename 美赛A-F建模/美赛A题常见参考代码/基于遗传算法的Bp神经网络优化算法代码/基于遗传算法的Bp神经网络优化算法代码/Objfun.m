function Obj=Objfun(X,P,T,hiddennum,P_test,T_test)
%% �����ֱ������Ⱥ�и��������Ŀ��ֵ
%% ����
% X�����и���ĳ�ʼȨֵ����ֵ
% P��ѵ����������
% T��ѵ���������
% hiddennum����������Ԫ��
% P_test:������������
% T_test:���������������
%% ���
% Obj�����и����Ԥ��������Ԥ�����ķ���

[M,N]=size(X);
Obj=zeros(M,1);
for i=1:M
%     Obj(i)=Bpfun(X(i,:),P,T,hiddennum,P_test,T_test);
    Obj(i)=BpFunction(X(i,:),P,T,hiddennum,P_test,T_test);
end