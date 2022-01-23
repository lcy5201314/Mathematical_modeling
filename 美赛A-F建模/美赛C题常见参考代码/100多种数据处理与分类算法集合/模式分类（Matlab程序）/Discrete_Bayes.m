function test_targets = Discrete_Bayes(train_patterns, train_targets, test_patterns, cost)

% Classify discrete patterns using the Bayes decision theory
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	cost			- Cost for class 0 (Optional, Unused yet)
%
% Outputs
%	test_targets	- Predicted targets

%First, find out if the patterns are indeed discrete 
%Of course, since this program is in Matlab, we define discrete as having 
%no more than one decimal place 
if (sum(sum(train_patterns*10~=floor(train_patterns*10))) ~= 0),
   error('Patterns are not discrete (See the definition of discrete in the m-file)')
end

Uc  = unique(train_targets);

%Find how the patterns are distributed
N	            = 0;
unique_patterns = [];
counts			= [];
for i = 1:size(train_patterns,2),
    data    = train_patterns(:,i);
    indices = find(sum((data*ones(1,size(train_patterns,2))-train_patterns).^2)==0);
    if isempty(unique_patterns),
        unique_patterns(:,1) = data;
        for j = 1:length(Uc),
            counts(j,1) = length(find(train_targets(indices) == Uc(j)));
        end
        N   = 2;
    else
        if isempty(find(sum((data*ones(1,size(unique_patterns,2))-unique_patterns).^2)==0)),
            %Add this pattern to the bank
            unique_patterns(:,N) = data;
            for j = 1:length(Uc),
                counts(j,N) = length(find(train_targets(indices) == Uc(j)));
            end
            N	= N + 1;
        end
    end
end

Px_given_w = (counts ./ (ones(length(Uc),1)* sum(counts)));
p          = hist(train_targets, Uc)' / length(train_targets);

Pw_given_x = Px_given_w .* (p*ones(1,N-1));
Pw_given_x = Pw_given_x ./ (ones(length(Uc),1)*sum(Pw_given_x));

[m, unique_targets] = max(Pw_given_x);
unique_targets      = unique_targets - 1;

test_targets = Nearest_Neighbor(unique_patterns,unique_targets,test_patterns,1);
