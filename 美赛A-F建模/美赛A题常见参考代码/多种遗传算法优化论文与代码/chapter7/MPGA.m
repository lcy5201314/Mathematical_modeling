%% ����Ⱥ�Ŵ��㷨
clear;
clc
close all
NIND=40;               %������Ŀ
NVAR=2;                %������ά��
PRECI=20;              %�����Ķ�����λ��
GGAP=0.9;             %����
MP=10;                   %��Ⱥ��Ŀ
FieldD=[rep(PRECI,[1,NVAR]);[-3,4.1;12.1,5.8];rep([1;0;1;1],[1,NVAR])];  %�������
for i=1:MP
    Chrom{i}=crtbp(NIND, NVAR*PRECI);                       %������ʼ��Ⱥ
end
pc=0.7+(0.9-0.7)*rand(MP,1);    %�ڡ�0.7,0.9����Χ����������������
pm=0.001+(0.05-0.001)*rand(MP,1);  %�ڡ�0.001,0.05����Χ����������������
gen=0;  %��ʼ�Ŵ�����
gen0=0; %��ʼ���ִ���
MAXGEN=10;  %���Ÿ������ٱ��ִ���
maxY=0; %����ֵ
for i=1:MP
    ObjV{i}=ObjectFunction(bs2rv(Chrom{i}, FieldD));%�������ʼ��Ⱥ�����Ŀ�꺯��ֵ
end
MaxObjV=zeros(MP,1);           %��¼������Ⱥ
MaxChrom=zeros(MP,PRECI*NVAR); %��¼������Ⱥ�ı���
while gen0<=MAXGEN
    gen=gen+1;       %�Ŵ�������1
    for i=1:MP
        FitnV{i}=ranking(-ObjV{i});                      % ����Ⱥ����Ӧ��
        SelCh{i}=select('sus', Chrom{i}, FitnV{i},GGAP); % ѡ�����
        SelCh{i}=recombin('xovsp',SelCh{i}, pc(i));      % �������
        SelCh{i}=mut(SelCh{i},pm(i));                    % �������
        ObjVSel=ObjectFunction(bs2rv(SelCh{i}, FieldD)); % �����Ӵ�Ŀ�꺯��ֵ
        [Chrom{i},ObjV{i}]=reins(Chrom{i},SelCh{i},1,1,ObjV{i},ObjVSel);    %�ز������
    end
    [Chrom,ObjV]=immigrant(Chrom,ObjV);     % �������
    [MaxObjV,MaxChrom]=EliteInduvidual(Chrom,ObjV,MaxObjV,MaxChrom);     % �˹�ѡ�񾫻���Ⱥ
    YY(gen)=max(MaxObjV);    %�ҳ�������Ⱥ�����ŵĸ���
    if YY(gen)>maxY   %�жϵ�ǰ�Ż�ֵ�Ƿ���ǰһ���Ż�ֵ��ͬ
        maxY=YY(gen); %��������ֵ
        gen0=0;
    else
        gen0=gen0+1; %����ֵ���ִ�����1
    end
end
%% ��������ͼ
plot(1:gen,YY)
xlabel('��������')
ylabel('���Ž�仯')
title('��������')
xlim([1,gen])
%% ������Ž�
[Y,I]=max(MaxObjV);    %�ҳ�������Ⱥ�����ŵĸ���
X=(bs2rv(MaxChrom(I,:), FieldD));   %���Ÿ���Ľ����
disp(['����ֵΪ��',num2str(Y)])
disp(['��Ӧ���Ա���ȡֵ��',num2str(X)])
