function test_targets = Local_Polynomial(train_patterns, train_targets, test_patterns, Nlp)

% Classify using the local polynomial fitting
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Nlp		        - Number of test points
%
% Outputs
%	test_targets	- Predicted targets

[M,L]	= size(train_patterns);
N       = size(test_patterns,2);

%Choose h
Ntest  = Nlp;
Ntrain = L - Ntest;
[train_indices, test_indices] = make_a_draw(Ntest, L);

h = 0;
for i = 1:Ntest,
    dist  = sum((train_patterns(:,train_indices) - train_patterns(:,test_indices(i))*ones(1,Ntrain)).^2);
    dist  = sort(dist);
    new_h = dist(round(Ntrain/10))/2;
    if (new_h > h),
        h  = new_h;
    end
end

%Classify all the test points to one of the Ntest points
Dn  = zeros(1, N);
for i = 1:N,
    dist        = sum((test_patterns(:,i) * ones(1,Ntest) - train_patterns(:,test_indices)).^2);
    [m, Dn(i)]  = min(dist);
end

%Now, built the plug-in classifier for each test point, and classify all the 
%points near it according to this classifier
test_targets = zeros(1, N);
for i = 1:Ntest,
   point   	      = train_patterns(:,test_indices(i));
   target_point   = train_targets(:,test_indices(i));
   theta   	   	  = fminunc('loglikelihood',zeros(M+1,1),optimset('Display','off'),train_patterns(:,train_indices),h,point,target_point);
   indices        = find(Dn == i);
   X		      = train_patterns(:,test_indices(i))*ones(1,length(indices))- test_patterns(:,indices);
   f_theta        = 1./(1+exp(-theta(1:end-1)'*X-theta(end)));
   test_targets(indices) = f_theta;
end

if (length(unique(train_targets)) == 2)
    test_targets = test_targets > 0.5;
end

