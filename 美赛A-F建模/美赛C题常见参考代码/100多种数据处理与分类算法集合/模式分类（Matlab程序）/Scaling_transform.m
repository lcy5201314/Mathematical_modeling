function [new_patterns, train_targets, var_mat, means] = Scaling_transform(train_patterns, train_targets, param, plot_on)

%Reshape the data points using the scaling transform
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	param			- Unused
%   plot_on         - Unused
%
%Outputs
%	new_patterns    - New patterns
%	targets			- New targets
%   var_mat			- Variance matrix
%   means           - Means vector

[r,c]		 = size(train_patterns);
means        = mean(train_patterns')';

new_patterns = train_patterns - means*ones(1,c);

var_mat      = var(new_patterns')';
new_patterns = new_patterns ./ (var_mat * ones(1, c));
