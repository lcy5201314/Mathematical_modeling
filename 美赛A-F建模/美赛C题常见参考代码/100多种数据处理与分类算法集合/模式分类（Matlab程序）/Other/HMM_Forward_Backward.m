function [a, b] = HMM_Forward_Backward(a, b, V)

% Find the probability transition matrices a,b from sample data using the forward-backward algorithm
%
% Inputs:
%	a					- Initial estimate of the transition probability matrix 
%	b					- Initial estimate of the output generator matrix 
%	V					- Observed output sequence
%
% Output:
%	a					- Estimated transition probability matrix 
%	b					- Estimated output generator matrix 

old_a = randn(size(a));
old_b = randn(size(b));

initial_alpha = (ones(1,size(a,1))/size(a,1));

theta = 1e-6;
while (max(sum(sum(abs(a-old_a))),sum(sum(abs(b-old_b)))) > theta),
    %Normalize a and b (for numeric reasons)
    a = a ./ (sum(a')' * ones(1,size(a,2)));
    b = b ./ (sum(b')' * ones(1,size(b,2)));

    %Compute alpha and beta
    [m, alpha] = HMM_forward (a,b,initial_alpha.*b(:,V(1))',V);
    [m, beta]  = HMM_backward(a,b,(ones(1,size(a,1))/size(a,1)).*b(:,V(end))',V);
    
    %Compute Gamma (via Zeta)
    for t = 1:length(V)-1,
        mul = alpha(:,t) * (beta(:,t+1).*b(:,V(t+1)))' .* a;
        if (sum(mul(:)) > 0),
            Zeta(:,:,t) = mul./sum(mul(:));
        end
        Gamma(:,t) = sum(Zeta(:,:,t),2);
    end
    Gamma(:,length(V)) = zeros(size(a,1),1);

    initial_alpha = Gamma(:,1)';
    
    %Compute a_hat from a and b by Eq. 140 (DHS Chapter 3)
    a_hat = sum(Zeta,3)./(sum(Gamma,2)*ones(1,size(a,2)));

    %Compute b_hat from a and b by Eq. 141 (DHS Chapter 3)
    b_hat = zeros(size(b));
    for j = 1:size(b,1),
        for k = 1:size(b,2),
            b_hat(j,k)=sum(Gamma(j,:).*(V==k))/sum(Gamma(j,:));
        end
    end
    
    %Save old a and b
    old_a = a;
    old_b = b;
    
    %a_ij <- a_hat_ij, b_jk <- b_hat_jk
    a = a_hat;
    b = b_hat;
end