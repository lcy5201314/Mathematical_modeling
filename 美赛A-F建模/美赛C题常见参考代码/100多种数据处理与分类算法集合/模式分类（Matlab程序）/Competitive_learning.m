function [patterns, targets, label, W] = Competitive_learning(train_patterns, train_targets, params, plot_on)

% Perform preprocessing using a competitive learning network
% Inputs:
% 	patterns	- Train patterns
%	targets	    - Train targets
%	params	    - [Number of partitions, learning rate]
%   plot_on     - Plot while performing processing?
%
% Outputs
%	patterns		- New patterns
%	targets			- New targets
%	label			- The labels given for each of the original patterns
%   W               - Weights matrice

max_iter       = 1000;
[c, r]		   = size(train_patterns);
[N, eta]	   = process_params(params);
decay          = 0.99;

%Preprocessing:
% x_i <- {x_i, 1}
x              = [train_patterns ; ones(1,r)];
%x_i <- x_i./||x_i||
x              = x ./ (ones(c+1,1) * sqrt(sum(x.^2)));

%Initialize the W's
i              = randperm(r);
W              = x(:,i(1:N));

for i = 1:max_iter,
    %Randomally order the patterns
    order = randperm(r);
    change= 0;
    
    for k = 1:r,
        J = W'*x(:,order(k));
        j = find(J == max(J));
        
        old_W   = W(:,j);
        
        %W_j <- W_j + eta*x
        W(:,j)  = W(:,j) + eta*x(:,order(k));
        
        %W_j <- W_j/||W_j||
        W(:,j)  = W(:,j) / sqrt(sum(W(:,j).^2));
        
        change = change + sum(abs(W(:,j) - old_W));
        
        if (plot_on > 0),
            %Assign each of the patterns to a center
            dist        = W'*x;
            [m, label]  = max(dist);
            centers     = zeros(c,N);
            for i = 1:N,
                in = find(label == i);
                if ~isempty(in)
                    centers(:,i) = mean(x(1:2,find(label==i))')';
                else
                    centers(:,i) = nan;
                end
            end
            
            %Plot the centers during the process 
            plot_process(centers, plot_on)
        end

    end

    eta = eta * decay;
    
    if (change/r < 1e-4),
        break
    end
    
end

if (i == max_iter),
   disp(['Maximum iteration (' num2str(max_iter) ') reached']);
else
    disp(['Finished after ' num2str(i) ' iterations.'])
end

%Assign each of the patterns to a center
dist        = W'*x;
[m, label]  = max(dist);
patterns     = zeros(c,N);
for i = 1:N,
    in = find(label == i);
    if ~isempty(in)
        patterns(:,i) = mean(x(1:end-1,find(label==i))')';
    else
        patterns(:,i) = nan;
    end
end

%Label the points
[m,label] = min(dist);
targets   = zeros(1,N);
Uc        = unique(train_targets);
for i = 1:N,
    n           = hist(train_targets(:,find(label == i)), Uc);
    [m, max_l]  = max(n);
    targets(i)  = Uc(max_l);
end