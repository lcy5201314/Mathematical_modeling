function targets = vote_C4_5(patterns, indices, tree, discrete_dim, Uc)  
    %Classify recursively using a tree  
    targets = zeros(1, size(patterns,2)); %设置每个样本的初始预测标签都是0 
      
    if (tree.dim == 0)  %这说明达到了树的叶子节点
        %Reached the end of the tree  
        targets(indices) = tree.child;  %得到样本对应的标签是tree.child
        return  
    end  
      
    %This is not the last level of the tree, so:  
    %First, find the dimension we are to work on  
    dim = tree.dim;  %得到分裂特征
    dims= 1:size(patterns,1);  %得到特征索引
      
    %And classify according to it  根据得到的决策树对测试样本进行分类
    if (discrete_dim(dim) == 0) %如果当前分裂特征是个连续特征 
        %Continuous pattern  
        in              = indices(find(patterns(dim, indices) <= tree.split_loc));  %找到当前测试样本中这个特征的特征值<=分裂值的样本索引
        targets     = targets + vote_C4_5(patterns(dims, :), in, tree.child(1), discrete_dim(dims), Uc);  %对这部分样本再分叉
        in              = indices(find(patterns(dim, indices) >  tree.split_loc));  %找到当前测试样本中这个特征的特征值>分裂值的样本索引
        targets     = targets + vote_C4_5(patterns(dims, :), in, tree.child(2), discrete_dim(dims), Uc);  %对这部分样本再分叉 
    else  %如果当前分裂特征是个离散特征
        %Discrete pattern  
        Uf              = unique(patterns(dim,:)); %得到这个样本集中这个特征的无重复特征值
        for i = 1:length(Uf)  %遍历每个特征值  
            if any(Uf(i) == tree.Nf)  %tree.Nf为树的分类特征向量 当前所有样本的这个特征的特征值
                in      = indices(find(patterns(dim, indices) == Uf(i)));  %找到当前测试样本中这个特征的特征值==分裂值的样本索引
                targets = targets + vote_C4_5(patterns(dims, :), in, tree.child(find(Uf(i)==tree.Nf)), discrete_dim(dims), Uc);%对这部分样本再分叉 
            end  
        end  
    end       
    %END 


