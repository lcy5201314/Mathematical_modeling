function [patterns, targets] = min_spanning_tree(train_patterns, train_targets, params, plot_on)

%Reduce the number of data points using a spanning tree
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	params			- [Method (NN-Nearest neighbor\inc - Inconsistant edge), Number of output data points or difference factor]
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets

if (nargin < 4),
    plot_on = 0;
end

[method, param]	= process_params(params);
[D,L]					= size(train_patterns);

%First, build a minimum spaning tree (MST)
%Stage 1: Compute the distance matrix
full_dist  = zeros(L);
for i = 1:L,
   full_dist(:,i) = sqrt(sum((train_patterns(:,i)*ones(1,L) - train_patterns).^2))';
end

%Stage 2: Connect the lines
dist	= zeros(L);
label = 1:L;

while (length(unique(label))>1),
   Ulabel		= unique(label);
   temp_dist	= zeros(length(Ulabel));   
   for i = 1:length(Ulabel)-1,
      for j = i+1:length(Ulabel),
         in_i = find(label == Ulabel(i));
         in_j = find(label == Ulabel(j));
         temp = full_dist(in_i, in_j);
         temp = temp(find(temp(:)>0));
         temp_dist(i, j) = min(min(temp));
      end
   end
   
   %Find which clusters to link
   ldist = sort(temp_dist(:));
   ldist = ldist(find(ldist~=0));
   [chosen_i, chosen_j] = find(temp_dist == ldist(1));
   
   %Now link them
   in_i 			= find(label == Ulabel(chosen_i));
   in_j 			= find(label == Ulabel(chosen_j));
   new_label	= min(label([in_i, in_j]));
   label(in_i)	= new_label;
   label(in_j) = new_label;
   [i, j]		= find(full_dist(in_i, in_j) == ldist(1));
   dist(in_i(i), in_j(j)) = full_dist(in_i(i), in_j(j));
   dist(in_j(j), in_i(i)) = full_dist(in_i(i), in_j(j));
end

%Now, process the MST according to the chosen method
if strcmp(method, 'NN'),
   %Nearest neighbor, so remove the longest edges in the graph   
   Dist = sort(dist(:));
   Dist = Dist(2:2:end);
   for i = 1:param,
      [j,k] = find(dist == Dist(end+1-i));
      j		= j(1);
      k		= k(1);
      dist(j,k) = 0;
      dist(k,j) = 0;
   end
      
else
   %Inconsistant edge
   temp_dist = dist;
   for i = 1:L,
      Dist	 	= dist(i,:);
      m_dist 	= mean(Dist(find(Dist~=0)));
      j			= find(Dist > param*m_dist);
      temp_dist(j,i)= 0;
   end
   dist = temp_dist;
end

%Now, make the MST into clusters 
label 		= zeros(1,L);
max_label	= 0;
while(~isempty(find(label==0))),
   in		= find(label == 0);  
   stack = in(1);
   
   if (label(stack) == 0),
	   max_label	= max_label + 1;
      curr_label = max_label;
   else
      curr_label = label(stack);
   end
         
   while (~isempty(stack)),
      %Label the first indice on the stack 
      label(stack(1)) = curr_label;
      
      %Add all those connected to this indice (and still unlabeled) to the stack
      j		= find(dist(stack(1),:) ~= 0);
		stack = [stack j(find(label(j)==0))];   
      
      %Remove this indice from the stack
      stack	= stack(2:end);
   end
end

%Build new patterns, and assign 0/1 labels to the original data
Ulabels  = unique(label);
Uc       = unique(train_targets);
Nclusters= length(Ulabels);
patterns = zeros(D, Nclusters);
targets  = zeros(1, Nclusters);
Labels	 = zeros(1, L);
for i = 1:Nclusters,
    in = find(label == Ulabels(i));
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
