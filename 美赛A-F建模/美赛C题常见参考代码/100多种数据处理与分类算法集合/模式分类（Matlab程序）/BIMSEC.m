function [patterns, targets, label, J] = BIMSEC(train_patterns, train_targets, params, plot_on)

%Reduce the number of data points using the basic iterative MSE clustering algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	params			- Algorithm parameters: [Number of output data points, Number of attempts]
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets
%	label			- The labels given for each of the original patterns

if (nargin < 4),
    plot_on = 0;
end
[Nmu, Ntries] = process_params(params);

[D,L]	= size(train_patterns);
dist	= zeros(Nmu,L);
label   = zeros(1,L);
Uc      = unique(train_targets);

%Initialize the mu's
mu			= randn(D,Nmu);
mu			= sqrtm(cov(train_patterns',1))*mu + mean(train_patterns')'*ones(1,Nmu);
ro          = zeros(1,Nmu);
n           = zeros(1,Nmu);
Ji          = zeros(1,Nmu);
J           = 1;
iter        = 1;

if (Nmu == 1),
   mu		= mean(train_patterns')';
   label	= ones(1,L);
else  
    %Assign each example to one of the mu's
    %Compute distances
    dist    = zeros(Nmu, L);
    for i = 1:Nmu,
        dist(i,:) = sqrt(sum((mu(:,i)*ones(1,L) - train_patterns).^2));
    end
    [m, label]  = min(dist);
    n           = hist(label, Nmu);

    while (Ntries > 0),
        
        iter        = iter + 1;
        J(iter)     = 0;
        
        %Select a sample x_hat  
        r     = randperm(L);
        x_hat = train_patterns(:,r(1));
        
        %i <- argmin||mi - x_hat||
        dist  = sqrt(sum((mu - x_hat * ones(1,Nmu)).^2));
        i     = find(dist == min(dist));
        
        %Compute ro if n(i) ~= 1
        if (n(i) ~=1),
            for j = 1:Nmu,
                if (i ~= j),
                    ro(j) = n(j)/(n(j)+1)*dist(j)^2;
                else
                    ro(j) = n(j)/(n(j)-1)*dist(j)^2;
                end
            end
            
            %Transfer x_hat if needed
            [m, k] = find(min(ro) == ro);
            if (k ~= i),
                label(r(1)) = k;
                n(i)        = n(i) - 1;
                n(k)        = n(k) + 1;
                
                %Recompute Je, and the mu's
                for j = 1:Nmu,
                    indexes = find(label == j);
                    mu(:,j) = mean(train_patterns(:,indexes)')';
                    Ji(j)   = sum(sum((mu(:,j)*ones(1,length(indexes)) - train_patterns(:,indexes)).^2));
                end
                
                J(iter)     = sum(Ji);
            end
                
        end
                 
        %Plot the centers during the process 
        plot_process(mu, plot_on)
 
        if (J(iter) == J(iter-1)),
            Ntries = Ntries - 1;
        end

    end
end

%Classify all the patterns to one of the mu's (1-NN)
dist = zeros(Nmu,L);
for i = 1:Nmu,
   dist(i,:) = sum((train_patterns - mu(:,i)*ones(1,L)).^2);
end
   
%Label the points
[m,label] = min(dist);
targets   = zeros(1,Nmu);
for i = 1:Nmu,
    N = hist(train_targets(:,find(label == i)), Uc);
    [m, max_l] = max(N);
    targets(i) = Uc(max_l);
end

patterns = mu;