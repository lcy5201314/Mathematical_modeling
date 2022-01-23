clear
clc

tic                                                                 %��ʼ��ʱ

[num_Citys,CityPosition]=ReadTSPFile('ulysses22.tsp');              %��ȡ.tsp�ļ�
%% ������������֮��ľ���
h=pdist(CityPosition);
D=squareform(h);
%% ��ʼ������
FishNum=9;                                                          %����10ֻ�˹���
Max_gen=200;                                                        %����������
trynumber=500;                                                      %�����̽����
Visual=16;                                                          %��֪����
deta=0.8;                                                           %ӵ��������
%% ��Ⱥ��ʼ��,ÿһ�б�ʾһ����
initFish=AF_init(FishNum,num_Citys);

BestX=zeros(Max_gen,num_Citys);                                     %��¼ÿ�ε�������������·��
BestY=zeros(Max_gen,1);                                             %��¼ÿ�ε�������������·���ľ���
besty=inf;                                                          %�����ܾ��룬��ʼ��Ϊ�����
gen=1;
currX=initFish;
currY=AF_foodconsistence(currX,D);
while gen<=Max_gen
    for i=1:FishNum
        [Xinext,flag]= AF_movestrategy(currX,i,D,Visual,deta,trynumber);
        currX(i,:)=Xinext;
    end
    currY=AF_foodconsistence(currX,D);
    [Ymin,index]=min(currY);
    if Ymin<besty
        besty=Ymin;
        bestx=currX(index,:);
        BestY(gen)=besty;
        BestX(gen,:)=bestx;
    else
        BestY(gen)=BestY(gen-1);
        BestX(gen,:)=BestX(gen-1,:);
    end
    disp(['��',num2str(gen),'�ε���,�ó�������ֵ��',num2str(BestY(gen))]);
    gen=gen+1;
    
end
figure
plot(1:Max_gen,BestY)
xlabel('��������')
ylabel('�Ż�ֵ')
title('��Ⱥ�㷨��������')
s=num2str(bestx(1));
for i=2:num_Citys
    s=strcat(s,'->');
    s=strcat(s,num2str(bestx(i)));
end
s=strcat(s,'->');
s=strcat(s,num2str(bestx(1)));
disp(['�ó�������·��:',s,',����ֵ:',num2str(besty)]);


toc                                                                 %������ʱ