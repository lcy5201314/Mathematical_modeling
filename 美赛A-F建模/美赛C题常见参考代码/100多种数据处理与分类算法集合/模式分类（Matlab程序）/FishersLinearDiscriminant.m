function [patterns, train_targets, w] = FishersLinearDiscriminant(train_patterns, train_targets, param, plot_on)

%Reshape the data points using the Fisher's linear discriminant
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	param			- Unused
%   plot_on         - Unused
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets
%   w				- Weights vector

train_one  = find(train_targets == 1);
train_zero = find(train_targets == 0);

s0			  = cov(train_patterns(:,train_zero)',1);
m0			  = mean(train_patterns(:,train_zero)');
s1			  = cov(train_patterns(:,train_one)',1);
m1			  = mean(train_patterns(:,train_one)');

sw			  = s0 + s1;
w			  = inv(sw)*(m0-m1)';
patterns      = [w'*train_patterns; zeros(1,length(train_targets))]; %We add a dimension because the toolbox needs 2D data
