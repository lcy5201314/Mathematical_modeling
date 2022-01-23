function [result] = statistics(tn, rnode, PValue, discrete_dim)  
    TypeName = {'1','2'};  
    TypeNum = [0 0];   
    test_patterns = PValue(:,1:end-1)';    
    class_num = length(TypeNum);
    type = zeros(tn,size(test_patterns,2));  
    for i = 1:tn  %对测试向量进行投票，共有tn棵树  
        type(tn,:) = vote_C4_5(test_patterns, 1:size(test_patterns,2), rnode{i,1}, discrete_dim, class_num); 
%         if strcmp( type(tn,:),TypeName{1})  
%             TypeNum(1) = TypeNum(1) + 1;  
%         else 
%             TypeNum(2) = TypeNum(2) + 1;   
%         end  
    end  
    result = mode(type,1)';%统计每列出现最多的数  分类问题
    
end  

 