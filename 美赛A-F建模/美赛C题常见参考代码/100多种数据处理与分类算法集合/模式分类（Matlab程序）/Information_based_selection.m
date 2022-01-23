function [patterns, targets, remaining_patterns] = Information_based_selection(patterns, targets, Npatterns)

% Koller and Sawami algorithm for pattern selection
%
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	Npatterns  		- Output dimension
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets
%	remaining_patterns	- The numbers of the selected features

Kdist			= 2;					%How many patterns to group together

%First, calculate the cross-entropy matrix
gamma			= infomat(patterns,targets);
Nf              = length(gamma);

gamma           = gamma + abs(min(min(gamma)));
diagonal		= diag(gamma);
gamma			= gamma - diag(diagonal);
discarded   = [];

%Discard redundant patterns
gamma       = gamma';
for k = 1:Nf-Npatterns,
   tic
   sgamma        = sort(gamma);
   sums			 = sum(sgamma(k:k+Kdist-1,:));
   sums(discarded) = inf;
   [m, min_i]	 = min(sums);
   discarded	 = [discarded min_i(1)];
   gamma(min_i(1),:) = 0;

   t = toc;
   disp([num2str(k) ':Discarded pattern number: ' num2str(min_i(1)) '. This took ' num2str(t) '[sec]'])
end

remaining_patterns = 1:Nf;
remaining_patterns(discarded) = 0;
remaining_patterns = remaining_patterns(find(remaining_patterns~=0));

disp(['Last two remaining pattern numbers: ' num2str(remaining_patterns)])

patterns = patterns(remaining_patterns, :);