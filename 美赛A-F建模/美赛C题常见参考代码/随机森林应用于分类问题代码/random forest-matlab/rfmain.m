clear all;  
rnode=cell(3,1);%3*1�ĵ�Ԫ����  
% rchild_value=cell(3,1);%3*1�ĵ�Ԫ����  
% rchild_node_num=cell(3,1);%3*1�ĵ�Ԫ����  
sn=300; %������ظ��ĳ�ȡsn������  
tn=10;  %ɭ���о���������Ŀ
load('aaa.mat');  
n = size(r,1);
%% ����ѵ���������ɭ�ֺ�ID3�㷨��������ɭ��  
discrete_dim = [];
for j=1:tn  
    Sample_num=randi([1,n],1,sn);%��1��1000�������ȡsn������  
    SData=r(Sample_num,:);  
    [tree,discrete_dim]= train_C4_5(SData, 5, 10, discrete_dim);  
    rnode{j,1}=tree;  
end  
      
%% ��������  
load('aaa.mat');  
T = r;
%TData=roundn(T,-1);  
TData = roundn(T,-1);  
%ͳ�ƺ�����������Ĳ�����������ͶƱ��Ȼ��ͳ�Ƴ�ѡƱ��ߵı�ǩ�������  
result = statistics(tn, rnode, TData, discrete_dim);  
gd = T(:,end);
len = length(gd);
count = sum(result==gd);
fprintf('����%d���������ж���ȷ����%d\n',len,count);
    
    
 
      
      
