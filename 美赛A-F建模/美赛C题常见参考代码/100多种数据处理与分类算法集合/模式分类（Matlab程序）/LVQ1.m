function [patterns, targets] = LVQ1(train_patterns, train_targets, Nmu, plot_on)

%Reduce the number of data points using linear vector quantization
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

alpha   = 0.9;
L		= length(train_targets);
dist	= zeros(Nmu,L);
label   = zeros(1,L);
Dim     = size(train_patterns, 1);

%Initialize the mu's
mu			= randn(Dim,Nmu);
mu			= sqrtm(cov(train_patterns',1))*mu + mean(train_patterns')'*ones(1,Nmu);
mu_target   = rand(1,Nmu)>0.5;
old_mu	    = zeros(Dim,Nmu);

while (sum(sum(abs(mu - old_mu))) > 0.1),
   old_mu = mu;
   
   %Classify all the patterns to one of the mu's
   for i = 1:Nmu,
      dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
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
         Q		  = ones(Dim,1) * (2*(train_targets(indices) == mu_target(i)) - 1);
         mu(:,i) = mu(:,i) + mean(((train_patterns(:,indices)-mu(:,i)*ones(1,length(indices))).*Q)')'*alpha;
      end
      
   end
   
   alpha = 0.95 * alpha;
   
   %Plot centers during training
   plot_process(mu, plot_on)

end

%Label the data
targets = zeros(1,Nmu);
Uc      = unique(train_targets);

for i = 1:Nmu,
    in				= find(label == i);
    if ~isempty(in),
        h            = hist(train_targets(in), Uc);
        [m, best]    = max(h);
        targets(i)	 = Uc(best);
        if length(in) == 1,
            patterns(:,i)	= train_patterns(:,in);
        else
            patterns(:,i)  = mean(train_patterns(:,in)')';
        end
    else
        patterns(:,i) = nan;
    end   
end
patterns = mu;
