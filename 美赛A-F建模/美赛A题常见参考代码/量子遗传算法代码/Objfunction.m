function [Y,X]=Objfunction(x,lenchrom)
%% Ŀ�꺯��
% ����     x�������Ʊ���
%   lenchrom���������Ķ�����λ��
% ���     Y��Ŀ��ֵ
%          X��ʮ������
bound=[-3.0 12.1;4.1 5.8];   % �����Ա����ķ�Χ
%% ��binary����ת����ʮ��������
X=bin2decFun(x,lenchrom,bound);
%% ������Ӧ��-����ֵ
Y=sin(4*pi*X(1))*X(1)+sin(20*pi*X(2))*X(2);
