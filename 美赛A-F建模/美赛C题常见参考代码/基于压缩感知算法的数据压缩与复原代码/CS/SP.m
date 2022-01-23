function [xr iter_num res] = SP(y, T, K,m)
% SAMP: Sparsity Adaptive Matching Pursuit algoritm for compressed sensing.
% For theoretical analysis, please refer to the paper : 
% Thong. T. Do, Lu Gan and Trac D. Tran ,"Sparsity Adaptive Matching
% Purusit for practical compressed sensing" available at http://dsp.ece.rice.edu/cs

% Written by Thong Do(thongdo@jhu.edu)
% Updated on July, 26th 2008

% parameter usage:
%   y: Mx1 observation vector 
%   T: MxN measurement matrix
%   K: any positive integer value not larger than sparsity
%   sigma: noise energy when sensing
%   xr: reconstructed sparse signal
%   iter_num: number of iterations

% Initialization
iter_num = 0; 
actset_size = K;
active_set = [];
res =y;

flag=1;
% while (norm(res)>5)
for i=1:m
    % candidate list
    [val idx] = sort(abs(T'*res), 'descend');        
    candidate_set = union(active_set, idx(1:actset_size));
    a=size(candidate_set,1);
%     if a>=2*K
%         break;
%     end
    % search candidate from x by ls
    hat_x=pinv(T(:,candidate_set))*y;
    [val idx] = sort(abs(hat_x), 'descend');
    new_active_set = candidate_set(idx(1:actset_size));
    new_res = y-T(:,new_active_set)*pinv(T(:,new_active_set))*y;
       
    % update residual and active set    
    res = new_res;
    active_set= new_active_set;
    
    iter_num = iter_num +1;
% if iter_num==5
%     break
% end
end 
iter_num=iter_num;
% reconstruction
N = size(T,2);
xr = zeros(N,1);
xr_active_set = pinv(T(:,active_set))*y;
xr(active_set) = xr_active_set;