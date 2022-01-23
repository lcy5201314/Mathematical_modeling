function [patterns, targets, pattern_numbers] = genetic_culling(patterns, targets, params)

% Culling type genetic algorithm for feature selection
%
% Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	params  		- Algorithm parameters: [%groups in each chromosome, Output dimension , classifier type, classifier params]
%
% Outputs:
%	patterns		- New patterns
%	targets			- New targets
%	pattern_numbers	- The numbers of the selected features

[Fremove, Ndims, Type, Params] = process_params(params);
Fremove   = 1/Fremove;

Ngenerations = 100;
Nfold        = 5;

if (floor(Ndims/2) ~= Ndims/2),
   error('Classifier output dimensions must be even.')
end

disp(['Using ' num2str(Ndims) ' dimensions per classifier'])
disp(['Fraction removed at each iteration is: 1/' num2str(Fremove)])

[Dims, Nf]  = size(patterns);

%Permute the examples
in          = randperm(Nf);
patterns    = patterns(:,in);
targets     = targets(in);

Ngroups     = floor(Dims/Ndims);
Nremove     = floor(Ngroups/Fremove);     %How many groups to remove at each iteration
L           = Ngroups * Ndims;
Pcorrect    = 0;

Lf          = floor(Nf/Nfold)*Nfold;
Fin         = reshape([1:Lf], Lf/Nfold, Nfold)';
train_indices= zeros(Nfold,Lf/Nfold*(Nfold-1));
test_indices= zeros(Nfold,Lf/Nfold);

for i=1:Nfold,
    train_indices(i,:)  = reshape(Fin([1:i-1,i+1:Nfold],:),1,Lf*(Nfold-1)/Nfold);    
    test_indices(i,:) = Fin(i,:);
end

%Make the initial, random, partition
indices     = reshape(randperm(L), Ndims, Ngroups);
correct     = inf*ones(Nfold, Ngroups);

%Start iterating
for i = 1:Ngenerations,
    
    disp(['Iteration ' num2str(i) ': Best classifier so far has ' num2str(Pcorrect*100) '%.']);
    
    %First, get Nfold cross validation of error (only for those groups that changed)
    for j = 1:Ngroups,
        if (~isfinite(correct(1,j))),
            for k = 1:Nfold,
                Ptargets        = feval(Type, squeeze(patterns(indices(:,j),train_indices(k,:))), targets(train_indices(k,:)), patterns(indices(:,j),test_indices(k,:)), Params);
    	        correct(k,j)    = mean(xor(Ptargets, targets(test_indices(k,:))));
            end
        end
        
    end
    
    %Find which Nremove groups had the worse performance, and which Nremove*2 had the best ones
    [m, in]     = sort(mean(correct));
    in_remove   = in(1:Nremove);
    in_temp		 = in(Nremove+1:Ngroups);
    in_temp		 = in_temp(randperm(length(in_temp)));
    
    %If there are not enough groups to choose from, add some using a random draw
    in_temp = ones(ceil(Nremove*2/length(in_temp)),1) * in_temp;
    in_temp	= in_temp(:);
    in_temp	= in_temp(randperm(length(in_temp)));
    in_best = in_temp(1:Nremove*2);
        

    if (max(m) > Pcorrect),
        best = indices(:,in(Ngroups));
        Pcorrect = max(m);
    end

    %Permute the best groups, and put them in place of the worse ones
    for j = 1:Nremove,
        g1  = in_best(j*2-1);
        g2  = in_best(j*2);
        gt  = in_remove(j);
        
        i1(randperm(Ndims))  = indices(:,g1);
        i2(randperm(Ndims))  = indices(:,g2);
        
        indices(1:Ndims/2,gt) = i1(1:Ndims/2)';
        indices(Ndims/2+1:Ndims,gt) = i2(1:Ndims/2)';
        
    end
    
    correct(:, in_remove) = inf.*ones(Nfold, Nremove);
end


%Assign the best pattern group to the output
patterns        = patterns(best,:);
pattern_numbers = best;

disp(['The best patterns are: ' num2str(best')])