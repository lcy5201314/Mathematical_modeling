function [patterns, targets, w] = DSLVQ(train_patterns, train_targets, Nmu, plot_on)

%Reduce the number of data points using distinction sensitive linear vector quantization 
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	Nmu				- Number of output data points
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets
%	w				- Weights vector

if (nargin < 4),
    plot_on = 0;
end

Ndim    = size(train_patterns, 1);
alpha   = 0.9;
beta	= 0.1;
L		= length(train_targets);
dist	= zeros(Nmu,L);
label   = zeros(1,L);

%Initialize the mu's
mu			= randn(Ndim,Nmu);
mu			= sqrtm(cov(train_patterns',1))*mu + mean(train_patterns')'*ones(1,Nmu);
mu_target   = rand(1,Nmu)>.5;
old_mu	    = zeros(Ndim,Nmu);

%Initialize the weight vector
w			= ones(size(train_patterns,1),1);

while (sum(sum(abs(mu - old_mu))) > 0.1),
   old_mu = mu;
   
   %Classify all the patterns to one of the mu's
   for i = 1:Nmu,
      dist(i,:) = sum(((w*ones(1,L)).*(train_patterns - mu(:,i)*ones(1,L))).^2);      
   end
      
   %For each sample, ...
   for i = 1:L,
      %Find the nearest neighbor classified correctly, and the nearest one classified
      %incorrectly
      d	= dist(:,i).*(mu_target'-.5)*2;
      dp = d;dn = d;
      dp(find(dp<0)) = nan;
      dn(find(dn>0)) = nan;
      ci = find(dp == min(dp));
      cj = find(dn == max(dn));
      if (isempty(ci) | isempty(cj)),
         break
      end
      di = abs(train_patterns(:,i) - mu(:,ci));
  	   dj = abs(train_patterns(:,i) - mu(:,cj));
      wn = (di-dj)/sum(abs(di-dj));
  	   nw	= w + beta*(wn - w);
     	nw(find(nw>1)) 	= 1;
      nw(find(nw<1e-4)) = 1e-4;      
      w	= nw./sum(abs(nw));
   end
      
   %Label the points
   [m,label] = min(dist);

   %Label the mu's
	for i = 1:Nmu,
   	if (length(train_targets(:,find(label == i))) > 0),
      	mu_target(i) = (sum(train_targets(:,find(label == i)))/length(train_targets(:,find(label == i))) > .5);
	   end
	end	
   
   %Recompute the mu's
   for i = 1:Nmu,
      indices = find(label == i);
      if ~isempty(indices),
         Q		  = ones(Ndim,1) * (2*(train_targets(indices) == mu_target(i)) - 1);
         mu(:,i) = mu(:,i) + mean(((train_patterns(:,indices)-mu(:,i)*ones(1,length(indices))).*Q)')'*alpha;
      end
      
   end
   
   alpha = 0.95 * alpha;
   beta	= 0.95 * beta;
   
   %Plot centers during training
   plot_process(mu, plot_on)

end

%Label the points
[m,label] = min(dist);
targets   = zeros(1,Nmu);
Uc        = unique(train_targets);
for i = 1:Nmu,
    N = hist(train_targets(:,find(label == i)), Uc);
    [m, max_l] = max(N);
    targets(i) = Uc(max_l);
end

patterns = mu;