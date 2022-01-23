function [patterns, targets, pattern_numbers] = Sequential_Feature_Selection(patterns, targets, params)

% Perform sequential feature selection
%
% Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	params  		- Algorithm parameters: [Forward/Backward, Output dimension , classifier type, classifier params]
%
% Outputs:
%	patterns		- New patterns
%	targets			- New targets
%	pattern_numbers	- The numbers of the selected features

[type, final_Ndims, Type, Params] = process_params(params);

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

%Generate initial partitions
switch type
case 'Forward'
    groups = combinations(1:Dims,2);
    st     = 2;
    step   = 1;
case 'Backward'
    groups = 1:Dims;
    st     = Dims;
    step   = -1;
otherwise
    error('Unknown type.')
end

%Start iterating
for i = st:step:final_Ndims,
    
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
    
    if (i == final_Ndims)
        break
    end
    
    %Generate new feature groups according to the selection method
    if (strcmp(type, 'Forward'))
        %Add a new feature
        all_features      = 1:Dims;
        possible_features = all_features(find(~ismember(all_features, best_features)));
        groups            = [possible_features', ones(length(possible_features),1)*best_features];
    else
        %Remove a feature
        groups            = combinations(1:length(best_features),length(best_features)-1);
        for j = 1:length(best_features),
            groups(j,:) = best_features(groups(j,:));
        end
    end
end

%Assign the best pattern group to the output
patterns        = patterns(best_features,:);
pattern_numbers = best_features;

disp(['The best patterns are: ' num2str(best_features)])