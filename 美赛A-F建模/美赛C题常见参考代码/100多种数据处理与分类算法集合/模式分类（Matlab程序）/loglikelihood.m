function ll = loglikelihood(theta, patterns, h, center_point, cp_target) 

% Used by the polynomial fitting algorithm

[c,r] = size(patterns);

patterns = center_point * ones(1,r) - patterns;

K = exp(-sum(patterns.^2)/(2*h))/(2*pi*h);
f = 1./(1 + exp(-theta(1:end-1)'*patterns - theta(end)));
if isempty(find(f==0))
   L  = log(f);
   ll = -sum(K.*(cp_target.*L+~cp_target.*(1-L)));
else
   ll = nan;
end
