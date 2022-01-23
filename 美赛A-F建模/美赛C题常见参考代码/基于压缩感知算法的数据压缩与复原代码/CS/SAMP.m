function [xr iter_num] =SAMP(y, Phi, step_size, sigma);
% SAMP: Sparsity Adaptive Matching Pursuit algoritm for compressed sensing.
% For theoretical analysis, please refer to the paper : 
% Thong. T. Do, Lu Gan and Trac D. Tran ,"Sparsity Adaptive Matching
% Purusit for practical compressed sensing" available at http://dsp.ece.rice.edu/cs

% Written by Thong Do(thongdo@jhu.edu)
% Updated on July, 26th 2008

% parameter usage:
%   y: Mx1 observation vector 
%   Phi: MxN measurement matrix
%   step_size: any positive integer value not larger than sparsity
%   sigma: noise energy when sensing
%   xr: reconstructed sparse signal
%   iter_num: number of iterations

% Initialization
iter_num = 0; 
actset_size = step_size;
active_set = [];
res = y;
stg_idx = 1; % stage index

while (norm(res)>sigma)

    % candidate list
    [val idx] = sort(abs(Phi'*res), 'descend');        
    candidate_set = union(active_set, idx(1:actset_size));

    % finalist
    [val idx] = sort(abs(pinv(Phi(:,candidate_set))*y), 'descend');
    new_active_set = candidate_set(idx(1:actset_size));
    new_res = y-Phi(:,new_active_set)*pinv(Phi(:,new_active_set))*y;

    if (norm(new_res) >= norm(res))   
        % shift into a new stage
        stg_idx = stg_idx + 1;
        actset_size = stg_idx*step_size;
       
    else          
        % update residual and active set    
        res = new_res;
        active_set= new_active_set;
        
    end 
    
    iter_num = iter_num +1;                        %whileµÄ´ÎÊý

end % loop

% reconstruction
N = size(Phi,2);
xr = zeros(N,1);
xr_active_set = pinv(Phi(:,active_set))*y;
xr(active_set) = xr_active_set;