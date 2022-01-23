function targets = vote_C4_5(patterns, indices, tree, discrete_dim, Uc)  
    %Classify recursively using a tree  
    targets = zeros(1, size(patterns,2)); %����ÿ�������ĳ�ʼԤ���ǩ����0 
      
    if (tree.dim == 0)  %��˵���ﵽ������Ҷ�ӽڵ�
        %Reached the end of the tree  
        targets(indices) = tree.child;  %�õ�������Ӧ�ı�ǩ��tree.child
        return  
    end  
      
    %This is not the last level of the tree, so:  
    %First, find the dimension we are to work on  
    dim = tree.dim;  %�õ���������
    dims= 1:size(patterns,1);  %�õ���������
      
    %And classify according to it  ���ݵõ��ľ������Բ����������з���
    if (discrete_dim(dim) == 0) %�����ǰ���������Ǹ��������� 
        %Continuous pattern  
        in              = indices(find(patterns(dim, indices) <= tree.split_loc));  %�ҵ���ǰ�����������������������ֵ<=����ֵ����������
        targets     = targets + vote_C4_5(patterns(dims, :), in, tree.child(1), discrete_dim(dims), Uc);  %���ⲿ�������ٷֲ�
        in              = indices(find(patterns(dim, indices) >  tree.split_loc));  %�ҵ���ǰ�����������������������ֵ>����ֵ����������
        targets     = targets + vote_C4_5(patterns(dims, :), in, tree.child(2), discrete_dim(dims), Uc);  %���ⲿ�������ٷֲ� 
    else  %�����ǰ���������Ǹ���ɢ����
        %Discrete pattern  
        Uf              = unique(patterns(dim,:)); %�õ������������������������ظ�����ֵ
        for i = 1:length(Uf)  %����ÿ������ֵ  
            if any(Uf(i) == tree.Nf)  %tree.NfΪ���ķ����������� ��ǰ�����������������������ֵ
                in      = indices(find(patterns(dim, indices) == Uf(i)));  %�ҵ���ǰ�����������������������ֵ==����ֵ����������
                targets = targets + vote_C4_5(patterns(dims, :), in, tree.child(find(Uf(i)==tree.Nf)), discrete_dim(dims), Uc);%���ⲿ�������ٷֲ� 
            end  
        end  
    end       
    %END 


