%Code written by Harel Avissar and Danny Bickson
%Suplamentary material for the paper: 
%Distributed fault identification via non-parametric belief propagation.
%By D. Bickson, H. Avissar, D. Dolev, S. P. Boyd, A. Ihler and D. Baron.
%Submitted for publication. Aug 2009.
%Code available from http://www.cs.huji.ac.il/labs/danss/p2p/gabp/
% This is the local opt function

function [Xamb Lamb rounds]=Round_and_Local(s,K,H,p,n,y,sigma,LOFlag)
    %% Rounding
    [x_sort,ind_x] = sort(s,'descend');
    x_cand = []; l_cand = [];
    %x_sort
    for i = 1:min(max(K, 3*n*p),n)
        x_cur = -ones(n,1);
        x_cur(ind_x(1:i)) = 1;
        l_cur = L(H, p, n, x_cur, y, sigma);
        x_cand = [x_cand x_cur]; l_cand = [l_cand l_cur]; 
    end
    [l_sort,ind_l] = sort(l_cand,'ascend');
    X_amb = x_cand(:,ind_l(1:K)) ;%get ambiguity set
    L_amb = l_sort(1:K);
    
    iter = 0;
    if LOFlag
        %%LocalOpt
        x_cand = X_amb;
       l_cand = L_amb;
        EXIT_FLAG = 0;
        while(~EXIT_FLAG)
            x_cur = X_amb(:,1);x_best = x_cur;
            iter = iter+1;
            for j = 1:n
                x_cur(j) = -x_cur(j);
                if ismember(x_cur', x_cand', 'rows')
                    x_cur(j) = -x_cur(j);
                else
                    x_cand = [x_cand x_cur];
                    l_cand = [l_cand L(H, p, n, x_cur, y, sigma)];
                    x_cur(j) = -x_cur(j);
                end
            end
            [l_sort,ind_l] = sort(l_cand,'ascend');
            X_amb = x_cand(:,ind_l(1:K)); %get ambiguity set
            L_amb = l_sort(1:K);
            if all(x_best==X_amb(:,1))
                EXIT_FLAG = 1;
            end
        end
    end
    Xamb = X_amb;
    Lamb = L_amb;
    rounds = iter;    
end
