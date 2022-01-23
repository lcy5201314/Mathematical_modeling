function [patterns, targets] = Deterministic_SA(train_patterns, train_targets, params, plot_on)

%Reduce the number of data points using the deterministic simulated annealing algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	params	    	- [Number of output data points, cooling rate (Between 0 and 1)]
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

T		= max(eig(cov(train_patterns',1)'))/2;    %Initial temperature
Tmin    = T/500;                                 %Stopping temperature

[d,L]	= size(train_patterns);
label   = zeros(1,L);
dist	= zeros(Nmu,L);
iter    = 0;
max_change = 1e-3;

%Init the inclusion matrix
inclusion_mat   = rand(Nmu, L);
inclusion_mat   = inclusion_mat ./ (ones(Nmu,1)*sum(inclusion_mat));

if (Nmu == 1),
    %Initialize the mu's
    mu			= mean(train_patterns')';
else
    %Initialize the P
    P   = rand(Nmu,L);
    P   = P ./ (ones(Nmu,1)*sum(P));
    
    while (T > Tmin),  
        iter    = iter + 1;
        T = epsi * T;
        
        for i = 1:L,
            %For each node (example):
            %Recompute the mu's
            for i = 1:Nmu,
               mu(:,i) = sum(((ones(d,1)*P(i,:)).*train_patterns)')'./(sum(P(i,:))); 
            end
            
            %Find the distances from mu's to patterns 
            for i = 1:Nmu,
               dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
            end
            dist = exp(-dist/T);
            %In this implementation, s_i is equal to dist!
               
            %Compute Gibbs distribution
            P = dist ./ (ones(Nmu,1) * sum(dist));
            if (~isfinite(sum(sum(P))))
                disp('P is infinite. Stopping.')
                break
            end
        end
        
        %Plot centers during training
        plot_process(mu, plot_on)
    
    end    
end


%Label the data
dist	= zeros(Nmu,L);
for i = 1:Nmu,
   dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
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


