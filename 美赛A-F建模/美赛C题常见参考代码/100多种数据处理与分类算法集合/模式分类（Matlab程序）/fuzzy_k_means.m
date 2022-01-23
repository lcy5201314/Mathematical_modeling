function [patterns, targets] = fuzzy_k_means(train_patterns, train_targets, Nmu, plot_on)

%Reduce the number of data points using the fuzzy k-means algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	Nmu				- Number of output data points
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets

if (nargin < 4),
    plot_on = 0;
end

b		= 2;
L		= length(train_targets);
dist	= zeros(Nmu,L);
label   = zeros(1,L);
Dim     = size(train_patterns,1);

%Initialize the mu's
mu		= randn(Dim,Nmu);
mu		= sqrtm(cov(train_patterns',1))*mu + mean(train_patterns')'*ones(1,Nmu);
old_mu	= zeros(Dim,Nmu);

%Initialize the P's
P		= randn(Nmu,L);
old_P	= zeros(Nmu,L);

while ((sum(sum(abs(mu - old_mu) > 1e-5)) + sum(sum(abs(P - old_P) > 1e-5)) > 0)),
   old_mu = mu;
   old_P  = P;
   
   %Classify all the patterns to one of the mu's
   for i = 1:Nmu,
      dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
   end
   
   %Recompute P's
   P = (1./dist).^(1/(b-1));
   P = P ./ (ones(Nmu,1) * sum(P));
   
   %Recompute the mu's
   P  = P.^b;
   mu = (train_patterns * P') ./ (ones(Dim,1)*sum(P'));

   %Plot centers during training
   plot_process(mu, plot_on)

end

%Classify the patterns
[m,label] = max(P);
targets   = zeros(1,Nmu);
Uc        = unique(train_targets);
for i = 1:Nmu,
    N = hist(train_targets(:,find(label == i)), Uc);
    [m, max_l] = max(N);
    targets(i) = Uc(max_l);
end

patterns = mu;