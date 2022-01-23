clear all;  
rnode=cell(3,1);%3*1的单元数组  
% rchild_value=cell(3,1);%3*1的单元数组  
% rchild_node_num=cell(3,1);%3*1的单元数组  
sn=300; %随机可重复的抽取sn个样本  
tn=10;  %森林中决策树的数目
load('aaa.mat');  
n = size(r,1);
%% 样本训练采用随机森林和ID3算法构建决策森林  
discrete_dim = [];
for j=1:tn  
    Sample_num=randi([1,n],1,sn);%从1至1000内随机抽取sn个样本  
    SData=r(Sample_num,:);  
    [tree,discrete_dim]= train_C4_5(SData, 5, 10, discrete_dim);  
    rnode{j,1}=tree;  
end  
      
%% 样本测试  
load('aaa.mat');  
T = r;
%TData=roundn(T,-1);  
TData = roundn(T,-1);  
%统计函数，对输入的测试向量进行投票，然后统计出选票最高的标签类型输出  
result = statistics(tn, rnode, TData, discrete_dim);  
gd = T(:,end);
len = length(gd);
count = sum(result==gd);
fprintf('共有%d个样本，判断正确的有%d\n',len,count);
    
    
 
      
      
