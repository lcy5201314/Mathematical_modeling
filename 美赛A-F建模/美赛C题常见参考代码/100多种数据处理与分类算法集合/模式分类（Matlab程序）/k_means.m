function [patterns, targets, label] = k_means(train_patterns, train_targets, Nmu, plot_on)

%Reduce the number of data points using the k-means algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	Nmu				- Number of output data points
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets
%	label			- The labels given for each of the original patterns

if (nargin < 4),
    plot_on = 0;
end

[D,L]	= size(train_patterns);
dist	= zeros(Nmu,L);
label   = zeros(1,L);

%Initialize the mu's
mu		= randn(D,Nmu);
mu		= sqrtm(cov(train_patterns',1))*mu + mean(train_patterns')'*ones(1,Nmu);
old_mu	= zeros(D,Nmu);

switch Nmu,
case 0,
    mu      = [];
    label   = [];
case 1,
   mu		= mean(train_patterns')';
   label	= ones(1,L);
otherwise
   while (sum(sum(abs(mu - old_mu) > 1e-5)) > 0),
      old_mu = mu;
      
      %Classify all the patterns to one of the mu's
      for i = 1:Nmu,
         dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
      end
      
      %Label the points
      [m,label] = min(dist);
      
      %Recompute the mu's
      for i = 1:Nmu,
         mu(:,i) = mean(train_patterns(:,find(label == i))')';
      end

      %Plot the centers during the process 
      plot_process(mu, plot_on)
   end
end
   
%Classify the patterns
targets   = zeros(1,Nmu);
Uc        = unique(train_targets);
for i = 1:Nmu,
    if (length(unique(train_targets(:,find(label == i)))) == 1)
        targets(i) = unique(train_targets(:,find(label == i)));
    else
        N = hist(train_targets(:,find(label == i)), Uc);
        if (~isempty(N))
            [m, max_l] = max(N);
            targets(i) = Uc(max_l);
        end
    end
end

patterns = mu;