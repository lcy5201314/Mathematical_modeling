function currentPathM=PM2VM(currentPathM,freq)
%% ��·�����������ٶȾ�����
rand_M=rand(length(currentPathM));
rand_M=(rand_M+rand_M')./2;
currentPathM_1=currentPathM.*freq;
currentPathM_2=currentPathM.*rand_M;
currentPathM=currentPathM_1>currentPathM_2;