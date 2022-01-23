%% ����
% x��һ������ĳ�ʼȨֵ����ֵ
% P��ѵ����������
% T��ѵ���������
% hiddennum����������Ԫ��
% P_test:������������
% T_test:���������������
%% ���
% err��Ԥ��������Ԥ�����ķ���
function [err,T_sim]=BpFunction(x,P,T,hiddennum,P_test,T_test)
inputnum=size(P,1);                             % �������Ԫ����
% hiddennum=2*inputnum+1;                           % ��������Ԫ����
outputnum=size(T,1);                                % �������Ԫ����
%% ���ݹ�һ��
[p_train,ps_train]=mapminmax(P,0,1);
p_test=mapminmax('apply',P_test,ps_train);
[t_train,ps_output]=mapminmax(T,0,1);
%% ��ʼ����BP����
net=newff(p_train,t_train,hiddennum);               %������Ϊhiddennum����Ԫ
%�趨�����������
net.trainParam.epochs=1000;
net.trainParam.goal=1e-3;
net.trainParam.lr=0.01;
net.trainParam.showwindow=false;                    %�߰�MATLABʹ�� ����ʾͼ�ο�
%% BP�������ʼȨֵ����ֵ
w1num=inputnum*hiddennum;                                           %����㵽�����Ȩֵ����
w2num=outputnum*hiddennum;                                          %�����㵽������Ȩֵ����
% x=2*rand(1,w1num+hiddennum+w2num+outputnum)-1;                      %�漴����Ȩֵ
W1=x(1:w1num);                                                      %��ʼ����㵽�������Ȩֵ
B1=x(w1num+1:w1num+hiddennum);                                      %������Ԫ��ֵ
W2=x(w1num+hiddennum+1:w1num+hiddennum+w2num);                      %�����㵽������Ȩֵ
B2=x(w1num+hiddennum+w2num+1:w1num+hiddennum+w2num+outputnum);      %�������ֵ
net.iw{1,1}=reshape(W1,hiddennum,inputnum);                         %Ϊ�����������㵽������Ȩֵ��ֵ
net.lw{2,1}=reshape(W2,outputnum,hiddennum);                        %Ϊ������������㵽�����Ȩֵ��ֵ
net.b{1}=reshape(B1,hiddennum,1);                                   %Ϊ�������������Ԫ��ֵ��ֵ
net.b{2}=reshape(B2,outputnum,1);                                   %Ϊ��������������ֵ��ֵ
%% ��ʼѵ��
net = train(net,p_train,t_train);
%% ��������
t_sim = sim(net,p_test);
T_sim = mapminmax('reverse',t_sim,ps_output);   %����һ��
err=norm(T_sim-T_test);                         %Ԥ��������Խ����ķ���������ԽС˵��Ԥ���Խ׼ȷ���������Ϊ0��˵��Ԥ�����ȫ׼ȷ

