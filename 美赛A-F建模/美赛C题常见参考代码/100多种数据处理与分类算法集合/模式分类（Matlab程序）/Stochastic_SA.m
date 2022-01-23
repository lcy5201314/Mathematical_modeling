function [patterns, targets] = Stochastic_SA(train_patterns, train_targets, params, plot_on)

%Reduce the number of data points using the stochastic simulated annealing algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	params	    	- [Number of output data points, cooling rate (Between 0 and 1)]
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns			- New patterns
%	targets			- New targets

if (nargin < 5),
    plot_on = 0;
end

%Parameters:
[Nmu, epsi]	 = process_params(params);
T		= max(eig(cov(train_patterns',1)'))/2;    %Initial temperature
Tmin    = 0.01;                                 %Stopping temperature

[d,L]	= size(train_patterns);
label   = zeros(1,L);
dist	= zeros(Nmu,L);
iter    = 0;
max_change = 1e-3;

%Initialize the mu's
mu			= mean(train_patterns')';

%Init the inclusion matrix
inclusion_mat   = rand(Nmu, L);
[m, i]          = max(inclusion_mat);
inclusion_mat   = zeros(Nmu,L);
for j = 1:L,
    inclusion_mat(i(j),j) = 1;
end

if (Nmu >= 1),
    while (T > Tmin),  
        iter    = iter + 1;
        index   = randperm(L);
        T = epsi * T;
        
        for i = 1:L,
            %Select a node (example) randomally. Poll all nodes once

            %Calculate the energy in this configuration: Ea <- 1/2*sum(w_ij*s_i*s_j)
            Ea = energy(train_patterns, inclusion_mat);
    
            %Change the configuration and see what the energy is
            config = inclusion_mat(:,index(i));
            change = rand(Nmu,1);
            [m, j] = max(~config.*change);
            new_inclusion_mat               = inclusion_mat;
            new_inclusion_mat(:,index(i))   = 0;
            new_inclusion_mat(j,index(i))   = 1;
            Eb     = energy(train_patterns, new_inclusion_mat);
            
            if (Eb < Ea),
                inclusion_mat = new_inclusion_mat;
            else
                if (exp(-(Eb-Ea)/T) > rand(1)),
                    inclusion_mat = new_inclusion_mat;
                end
            end
        end
        
        %Recalculate the mu's
        mu  = zeros(d, Nmu);
        for i = 1:Nmu,
            indices = find(inclusion_mat(i,:) == 1);
            mu(:,i) = mean(train_patterns(:,indices)')';
        end

        %Plot centers during training
        plot_process(mu, plot_on)
    
    end    
end


%Classify new centers
dist	= zeros(Nmu,L);
for i = 1:Nmu,
   dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
end
[m,label] = min(dist);

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



function E = energy(patterns, inclusion_matrix)

%Calculate the energy value given the patterns and the inclusion matrix
%The energy function tries to minimize the in-class variance

[N,M]   = size(inclusion_matrix);
e       = zeros(1,N+1);

for i = 1:N,
    indices = find(inclusion_matrix(i,:) == 1);
    mu      = mean(patterns(:,indices)')';
    e(i)    = sum(sum((patterns(:,indices) - mu*ones(1,length(indices))).^2));
end

E = sum(e);