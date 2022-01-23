function test_targets = CART(train_patterns, train_targets, test_patterns, params)

% Classify using classification and regression trees
% Inputs:
% 	training_patterns   - Train patterns
%	training_targets	- Train targets
%   test_patterns       - Test  patterns
%	params		        - [Impurity type, Percentage of incorrectly assigned samples at a node]
%                         Impurity can be: Entropy, Variance (or Gini), or Misclassification
%
% Outputs
%	test_targets        - Predicted targets


[Ni, M]		   = size(train_patterns);

%Get parameters
[split_type, inc_node] = process_params(params);

%Build the tree recursively
disp('Building tree')
tree            = make_tree(train_patterns, train_targets, M, split_type, inc_node);

%Classify test patterns according to the tree
disp('Classify test patterns using the tree')
test_targets    = use_tree(test_patterns, 1:size(test_patterns,2), tree);

%END

function targets = use_tree(patterns, indices, tree)
%Classify recursively using a tree

if isnumeric(tree.Raction)
   %Reached an end node
   targets = zeros(1,size(patterns,2));
   targets(indices) = tree.Raction(1);
else
   %Reached a branching, so:
   %Find who goes where
   in_right    = indices(find(eval(tree.Raction)));
   in_left     = indices(find(eval(tree.Laction)));
   
   Ltargets	 = use_tree(patterns, in_left, tree.left);
   Rtargets	 = use_tree(patterns, in_right, tree.right);
   
   targets		 = Ltargets + Rtargets;
end
%END use_tree 

function tree = make_tree(patterns, targets, Dlength, split_type, inc_node)
%Build a tree recursively

if (length(unique(targets)) == 1),
   %There is only one type of targets, and this generates a warning, so deal with it separately
   tree.right      = [];
   tree.left       = [];
   tree.Raction    = targets(1);
   tree.Laction    = targets(1);
   return
end

[Ni, M] = size(patterns);
Nt      = unique(targets);
N       = hist(targets, Nt);

if ((sum(N < Dlength*inc_node) == length(Nt) - 1) | (M == 1)), 
   %No further splitting is neccessary
   tree.right      = [];
   tree.left       = [];
   if (length(Nt) ~= 1),
      MLlabel			  = find(N == max(N));
   else
      MLlabel 		  = 1;
   end
   tree.Raction    = Nt(MLlabel);
   tree.Laction    = Nt(MLlabel);
   
else
   %Split the node according to the splitting criterion  
   deltaI		= zeros(1,Ni);
   split_point = zeros(1,Ni);
   op				= optimset('Display', 'off');   
   for i = 1:Ni,
      split_point(i) = fminbnd('CARTfunctions', min(patterns(i,:)), max(patterns(i,:)), op, patterns, targets, i, split_type);
      I(i)				= feval('CARTfunctions', split_point(i), patterns, targets, i, split_type);
   end
   
   [m, dim]				= min(I);
   loc					= split_point(dim);
    
   %So, the split is to be on dimention 'dim' at location 'loc'
   indices		 = 1:M;
   tree.Raction= ['patterns(' num2str(dim) ',indices) >  ' num2str(loc)];
   tree.Laction= ['patterns(' num2str(dim) ',indices) <= ' num2str(loc)];
   in_right    = find(eval(tree.Raction));
   in_left     = find(eval(tree.Laction));
   
   if isempty(in_right) | isempty(in_left)
      %No possible split found
	   tree.right      = [];
   	tree.left       = [];
	   if (length(Nt) ~= 1),
   	   MLlabel  	 = find(N == max(N));
	   else
   	   MLlabel 		 = 1;
	   end
   	tree.Raction    = Nt(MLlabel);
   	tree.Laction    = Nt(MLlabel);
   else
	   %...It's possible to build new nodes
   	tree.right = make_tree(patterns(:,in_right), targets(in_right), Dlength, split_type, inc_node);
   	tree.left  = make_tree(patterns(:,in_left), targets(in_left), Dlength, split_type, inc_node);      
   end
   
end

