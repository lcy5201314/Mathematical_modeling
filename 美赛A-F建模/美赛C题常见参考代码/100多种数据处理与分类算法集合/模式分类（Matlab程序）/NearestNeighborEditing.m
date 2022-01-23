function [patterns, targets] = NearestNeighborEditing(train_patterns, train_targets, Nmu, plot_on)

%Reduce the number of data points using the nearest neighbor editing algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	Nmu				- Unused
%   plot_on         - Unused
%
%Outputs
%	patterns		- New patterns (Prototypes)
%	targets			- New targets

N   = size(train_patterns,2);
options = optimset('Display','off');

% In this implementation we find only adjacent Voronoi cells to each example, in order to reduce 
% the computational costs. This is done through linear programming (LP)
b   = sum(train_patterns.^2);
A   = 2*train_patterns;
A(end+1, :) = -ones(1,N);

neighbors = zeros(N);
for i = 1:N-1,
    for j = i+1:N,
        b_tag       = b; 
        b_tag(j)    = b(j) + 1;
        constraintA = [A, -A(:,i)]';
        constraintB = [b_tag, -b(i)]';
        sol         = linprog(-A(:,j)', constraintA, constraintB, [], [], [], [], [], options);
        f           = b(j) - A(:,j)'*sol;
        
        if (f < 0)
            neighbors(i,j) = 1;
            neighbors(j,i) = 1;
        end
    end
    disp (['Finished ' num2str(i) ' of ' num2str(N-1) ' loops.'])
end
        
% For each pattern, check if all it's neighbors are the same. If so, delete it.
keep = ones(1,N);
for i = 1:N,
    in  = find(neighbors(i,:) == 1);
    if (length(unique(train_targets(in))) == 1)
        keep(i) = 0;
    end
end

keep_in = find(keep);

patterns = train_patterns(:,keep_in);
targets  = train_targets (keep_in);
