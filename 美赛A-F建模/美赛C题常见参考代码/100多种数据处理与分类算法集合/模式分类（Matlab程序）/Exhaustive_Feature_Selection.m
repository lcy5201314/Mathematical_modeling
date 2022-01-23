function [patterns, targets, pattern_numbers] = Exhaustive_Feature_Selection(patterns, targets, params)

% Perform exhaustive (Brute-force) feature selection
%
% Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	params  		- Algorithm parameters: [Output dimension , classifier type, classifier params]
%
% Outputs:
%	patterns		- New patterns
%	targets			- New targets
%	pattern_numbers	- The numbers of the selected features

[final_Ndims, Type, Params] = process_params(params);

Nfold        = 5;
[Dims, Nf]   = size(patterns);

if (final_Ndims < 2),
   error('The minimum feature number is two.')
end

if (final_Ndims > Dims)
    error('The output dimension is larger than the input dimension.');
end

Lf           = floor(Nf/Nfold)*Nfold;
Fin          = reshape([1:Lf], Lf/Nfold, Nfold)';
train_indices= zeros(Nfold,Lf/Nfold*(Nfold-1));
test_indices = zeros(Nfold,Lf/Nfold);
for i=1:Nfold,
    train_indices(i,:)  = reshape(Fin([1:i-1,i+1:Nfold],:),1,Lf*(Nfold-1)/Nfold);    
    test_indices(i,:) = Fin(i,:);
end

%Generate all possible feature groups 
groups = combinations(1:Dims, final_Ndims);

%First, get Nfold cross validation of error 
Ng    = size(groups, 1);
score = zeros(Nfold, Ng);
for j = 1:Ng,
    for k = 1:Nfold,
        Ptargets        = feval(Type, patterns(groups(j,:),train_indices(k,:)), targets(train_indices(k,:)), patterns(groups(j,:),test_indices(k,:)), Params);
        score(k,j)      = 1 - mean(xor(Ptargets, targets(test_indices(k,:))));
    end
end

%Find which feature group gave the best results
ave_score = mean(score);
[max_score, best_group] = max(ave_score);
best_features = sort(groups(best_group,:));

%Assign the best pattern group to the output
patterns        = patterns(best_features,:);
pattern_numbers = best_features;

disp(['The best patterns are: ' num2str(best_features) ' (With ' num2str(100*max_score) '% success)'])