function [patterns, targets] = LVQ3(train_patterns, train_targets, Nmu, plot_on)

%Reduce the number of data points using linear vector quantization
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	Nmu				- Number of output data points
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns			- New patterns
%	targets			- New targets

if (nargin < 4),
    plot_on = 0;
end

if (length(unique(train_targets)) < 2),
    error('LVQ3 works only if there are patterns from at least two classes.')
end

L		= length(train_targets);
alpha = 10/L;
[D,r] = size(train_patterns);
dist	= zeros(Nmu,L);
label = zeros(1,L);

window = 0.25;
epsilon= 0.25;

%Initialize the mu's
mu			= randn(D,Nmu);
mu			= sqrtm(cov(train_patterns',1))*mu + mean(train_patterns')'*ones(1,Nmu);
mu_target= [zeros(1,floor(Nmu/2)) ones(1,Nmu-floor(Nmu/2))];
old_mu	= zeros(D,Nmu);

iterations = 0;

while ((sum(sum(abs(mu - old_mu))) > 0.01) & (iterations < 1e4)),
   iterations = iterations + 1;
   old_mu = mu;
   
   %Classify all the patterns to one of the mu's
   for i = 1:Nmu,
      dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
   end
   
   %Label the points
   [dist,label] = sort(dist);
   closest		 = dist(1:2,:);
    
   %Compute windows
   in_window = (min(closest(1,:)./closest(2,:), closest(2,:)./closest(1,:)) > (1-window)/(1+window));
   indices	 = find(in_window);
   
   %Move the mu's
   for i = 1:length(indices),
      x	 = indices(i);
      mu1 = label(1,x);
      mu2 = label(2,x);
      if ((train_targets(x) == mu_target(mu1)) & (train_targets(x) == mu_target(mu2))),
         mu(:,mu1) = mu(:,mu1) + epsilon * alpha * (train_patterns(:,x) - mu(:,mu1));
         mu(:,mu2) = mu(:,mu2) + epsilon * alpha * (train_patterns(:,x) - mu(:,mu2));
      else
         if (train_targets(x) == mu_target(mu1)),
	         mu(:,mu1) = mu(:,mu1) + alpha * (train_patterns(:,x) - mu(:,mu1));
            mu(:,mu2) = mu(:,mu2) - alpha * (train_patterns(:,x) - mu(:,mu2));
         else
	         mu(:,mu1) = mu(:,mu1) - alpha * (train_patterns(:,x) - mu(:,mu1));
            mu(:,mu2) = mu(:,mu2) + alpha * (train_patterns(:,x) - mu(:,mu2));
         end
      end
   end

   alpha = 0.95 * alpha;

   %Plot centers during training
   plot_process(mu, plot_on)

end

%Label the data
targets  = mu_target;
patterns = mu;

