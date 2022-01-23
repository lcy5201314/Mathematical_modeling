function [patterns, targets] = deterministic_annealing(train_patterns, train_targets, params, plot_on)

%Reduce the number of data points using the deterministic annealing algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	params			- [Number of output data points, Cooling rate]
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets

if (nargin < 4),
    plot_on = 0;
end

%Parameters:
[Nmu, epsi] = process_params(params);
T		    = max(eig(cov(train_patterns',1)'))/2;    %Initial temperature
Tmin        = 0.01;                                 %Stopping temperature

[d,L]	    = size(train_patterns);
Ncent       = 1;
label       = zeros(1,L);
dist	    = zeros(Ncent,L);
iter        = 0;
max_change  = 1e-3;

%Initialize the mu's
mu			= mean(train_patterns')';

while (T > Tmin),  
    iter = iter + 1;
    
    %Find the distances from mu's to patterns 
    for i = 1:Ncent,
       dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
    end
    dist = exp(-dist/T);
   
    %Compute Gibbs distribution
    P = dist ./ (ones(Ncent,1) * sum(dist));
   
    %Recompute the mu's
    old_mu = mu;   
    for i = 1:Ncent,
       mu(:,i) = sum(((ones(d,1)*P(i,:)).*train_patterns)')'./(sum(P(i,:))); 
    end

    if (sum(sum(abs(old_mu-mu))) <= max_change)
        %Minimum reached, so decrease temperature ...
        T = epsi * T;
        if (Ncent >= Nmu),
            %There are enough partitions
            break
        end
        
        %...and add a center near the center that has the most variance
        if (Ncent > 1),
            for i = 1:Ncent,
               dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
            end
            [m,label] = min(dist);
            %Find the variance of the patterns around all the centers
            Smu = zeros(1,Nmu);
            for i = 1:Ncent,
                Smu(i) = sum(std(train_patterns(:,find(label == i))'));
            end
            [m, max_std]  = max(Smu);
        else
            max_std = 1;
        end
        Ncent         = Ncent + 1;
        mu(:,Ncent)   = mu(:,max_std) + randn(d,1).*std(mu')'/10;
        mu(:,max_std) = mu(:,max_std) + randn(d,1).*std(mu')'/10;
        %mu(:,Ncent) = randn(2,1);
    end
    
    %Plot centers during training
    plot_process(mu, plot_on)
end

%Label the data
dist	= zeros(Ncent,L);
for i = 1:Ncent,
   dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
end
[m,label] = min(dist);

%Label the points
[m,label] = min(dist);
targets   = zeros(1,Ncent);
Uc        = unique(train_targets);
for i = 1:Ncent,
    N = hist(train_targets(:,find(label == i)), Uc);
    [m, max_l] = max(N);
    targets(i) = Uc(max_l);
end

patterns = mu;

