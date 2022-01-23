function test_targets = ID3(train_patterns, train_targets, test_patterns, params)

% Classify using Quinlan's ID3 algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	params		    - [Number of bins for the data, Percentage of incorrectly assigned samples at a node]
%
% Outputs
%	test_targets	- Predicted targets

[Ni, M]		   = size(train_patterns);

%Get parameters
[Nbins, inc_node] = process_params(params);
inc_node          = inc_node*M/100;

%First, bin the data and the decision region data
[H, binned_train_patterns, region]  = high_histogram(train_patterns, Nbins);
[H, binned_test_patterns]           = high_histogram(test_patterns, Nbins, region);

%Build the tree recursively
disp('Building tree')
tree        = make_tree(binned_train_patterns, train_targets, inc_node, Nbins);

%Classifying test patterns according to the tree
disp('Classifying test patterns using the tree')
test_targets    = use_tree(binned_test_patterns, 1:size(test_patterns,2), tree, Nbins, unique(train_targets));

%END

function targets = use_tree(patterns, indices, tree, Nbins, Uc)
%Classify recursively using a tree

targets = zeros(1, size(patterns,2));

if (isempty(indices))
    return;
end

if ((size(patterns,1) == 1) | (tree.split_dim == 0)),
    %Only one dimension left, so work on it
    for i = 1:Nbins,
        in = indices(find(patterns(indices) == i));
        if ~isempty(in),
            if isfinite(tree.child(i)),
                targets(in) = tree.child(i);
            else
                %No data was found in the training set for this bin, so choose it randomally
                n           = 1 + floor(rand(1)*length(Uc));
                targets(in) = Uc(n);
            end
        end
    end
    return
end
        
%This is not the last level of the tree, so:
%First, find the dimension we are to work on
dim = tree.split_dim;
dims= find(~ismember(1:size(patterns,1), dim));

%And classify according to it
for i = 1:Nbins,
    in      = indices(find(patterns(dim, indices) == i));
    targets = targets + use_tree(patterns(dims, :), in, tree.child(i), Nbins, Uc);
end
    
%END use_tree 

function tree = make_tree(patterns, targets, inc_node, Nbins)
%Build a tree recursively

[Ni, L]     = size(patterns);
Uc          = unique(targets);

%When to stop: If the dimension is one or the number of examples is small
if ((Ni == 1) | (inc_node > L)),
    %Compute the children non-recursively
    for i = 1:Nbins,
        tree.split_dim  = 0;
        if ~isempty(targets)
            indices         = find(targets == i);
            if ~isempty(indices),
                if (length(unique(targets(indices))) == 1),
                    tree.child(i) = targets(indices(1));
                else
                    H               = hist(targets(indices), Uc);
                    [m, T]          = max(H);
                    tree.child(i)   = Uc(T);
                end
            else
                tree.child(i)   = inf;
            end
        else
            tree.child(i)   = inf;
        end
    end
    return
end

%Compute the node's I
for i = 1:length(Uc),
    Pnode(i) = length(find(targets == Uc(i))) / L;
end
Inode = -sum(Pnode.*log(Pnode)/log(2));

%For each dimension, compute the gain ratio impurity
delta_Ib    = zeros(1, Ni);
P           = zeros(length(Uc), Nbins);
for i = 1:Ni,
    for j = 1:length(Uc),
        for k = 1:Nbins,
            indices = find((targets == Uc(j)) & (patterns(i,:) == k));
            P(j,k)  = length(indices);
        end
    end
    Pk          = sum(P);
    P           = P/L;
    Pk          = Pk/sum(Pk);
    info        = sum(-P.*log(eps+P)/log(2));
    delta_Ib(i) = (Inode-sum(Pk.*info))/-sum(Pk.*log(eps+Pk)/log(2));
end

%Find the dimension minimizing delta_Ib 
[m, dim] = max(delta_Ib);

%Split along the 'dim' dimension
tree.split_dim = dim;
dims           = find(~ismember(1:Ni, dim));
for i = 1:Nbins,
    indices       = find(patterns(dim, :) == i);
    tree.child(i) = make_tree(patterns(dims, indices), targets(indices), inc_node, Nbins);
end




