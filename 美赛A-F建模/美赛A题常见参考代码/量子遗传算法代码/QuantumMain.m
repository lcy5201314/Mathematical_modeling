clc;
clear all;
close all;
%----------------��������-----------------------
MAXGEN=200;                        % ����Ŵ�����
sizepop=40;                        % ��Ⱥ��С
lenchrom=[20 20];          % ÿ�������Ķ����Ƴ���
trace=zeros(1,MAXGEN);
%--------------------------------------------------------------------------
% ��Ѹ��� ��¼����Ӧ��ֵ��ʮ����ֵ�������Ʊ��롢���ӱ��ر���
best=struct('fitness',0,'X',[],'binary',[],'chrom',[]);   
%% ��ʼ����Ⱥ
chrom=InitPop(sizepop*2,sum(lenchrom));
%% ����Ⱥʵʩһ�β��� �õ������Ʊ���
binary=collapse(chrom); 
%% ����Ⱥ�������Ӧ��ֵ���Ͷ�Ӧ��ʮ����ֵ
[fitness,X]=FitnessFunction(binary,lenchrom);   % ʹ��Ŀ�꺯��������Ӧ��
%% ��¼��Ѹ��嵽best
[best.fitness bestindex]=max(fitness);     % �ҳ����ֵ
best.binary=binary(bestindex,:);
best.chrom=chrom([2*bestindex-1:2*bestindex],:);
best.X=X(bestindex,:);
trace(1)=best.fitness;
fprintf('%d\n',1)
%% ����
for gen=2:MAXGEN
    fprintf('%d\n',gen)  %��ʾ��������
    %% ����Ⱥʵʩһ�β���
    binary=collapse(chrom);
    %% ������Ӧ��
    [fitness,X]=FitnessFunction(binary,lenchrom);
    %% ������ת��
    chrom=Qgate(chrom,fitness,best,binary);
    [newbestfitness,newbestindex]=max(fitness);    % �ҵ����ֵ
    % ��¼��Ѹ��嵽best
    if newbestfitness>best.fitness
        best.fitness=newbestfitness;
        best.binary=binary(newbestindex,:);
        best.chrom=chrom([2*newbestindex-1:2*newbestindex],:);
        best.X=X(newbestindex,:);
    end
    trace(gen)=best.fitness;
end

%% ����������
plot(1:MAXGEN,trace);
title('��������');
xlabel('��������');
ylabel('ÿ���������Ӧ��');

%% ��ʾ�Ż����
disp(['���Ž�X��',num2str(best.X)])
disp(['���ֵY:',num2str(best.fitness)]);
