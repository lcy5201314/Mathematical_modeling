function [test_targets, best_fun] = genetic_programming(train_patterns, train_targets, test_patterns, params)

% A genetic programming algorithm for classification
%
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	params  		- [Initial function length, Number of generations, Number of solutions]
%
%Outputs
%	test_targets	- Predicted targets
%   func			- The best solution found

IterDisp  = 10;

[Ni, M]   = size(train_patterns);
x		  = train_patterns;
Uc        = unique(train_targets);

%Find parameters
[Flength, Ngenerations, Nsol] = process_params(params);

if (Flength) > 20,
   warning('This may take some time!')
end

%Define basic operators. If division by zero is unwanted, remove the devision operator
warning off %To suppress all those divide by zero warnings
operators	 = [' +'; ' -'; '.*'; './'];
base			 = size(operators,1);
for i = 1:Ni,
   operators(base+i,:) = ['X' num2str(i)];
end

%Define the basic tree shape
leafs = [];
solutions = zeros(1,Nsol); %Pointer to the first leaf in the solution

for i = 1:Nsol;
   solutions(i) = length(leafs) + 1;
   leafs			 = foliage(leafs, operators, Flength);
end

disp(['Starting with ' num2str(Nsol) ' solutions, totaling ' num2str(length(leafs)) ' leaves.'])

best = 0;

%Start doing the mutations
for i = 1:Ngenerations,
   %Rank each solution
   ranking = zeros(1,Nsol);
   for j = 1:Nsol,
      func    = collect_tree(solutions(j), leafs);
      for k = 1:Ni,
         func = strrep(func, ['X' num2str(k)], ['x(' num2str(k) ',:)']);
      end
      
      %Check if this string is a valid solution, and if so, put the ranking into place
      try
         t = eval(func);
         if (length(Uc) == 2)
             t = t > 0;
         else
             t = round(t);
         end
         ranking(j) = sum(t == train_targets) / length(train_targets);
      catch
      end
      
   end
   
   %Display
   if (i/IterDisp == floor(i/IterDisp)),
      disp(['Iteration ' num2str(i) ': Best solution has ' num2str(max(ranking)*100) '% correct (The best so far was ' num2str(best*100) '% correct).'])
   end
      
   %Sort the solutions according to rank
   [m, indices] = sort(ranking);
   indices		 = fliplr(indices);
   
   %It is sometimes good to save the best solution so far, in case the mutations ruin it!
   if (max(ranking) > best),
      best 		= max(ranking);
      best_fun = collect_tree(solutions(indices(1)), leafs);
   end
   
   %Start the genetic operations
   for j = 1:floor(Nsol/2),
      %Select two solutions according to their rank
      s1	= indices(j*2-1);
      s2 = indices(j*2);
      
      %Now do the genetic operators
      %There is no insertion in this implementation, as it generates too many impossible solutions
      
      k	= 1 + floor(3*rand(1));
      switch k,
      case 1,
         %Replications
         %Replicate from 2 to 1 or vice versa?
         if rand(1) > .5,
            temp = s1;
            s1   = s2;
            s2   = temp;
         end
         
         %Find where from to replicate
         node_number1 = walk_along(solutions(s1), leafs);
         node_number2 = walk_along(solutions(s2), leafs);
         [new_node, leafs] = copy(node_number1, leafs);
         
         %Put on left or right?
         if rand(1) > .5,
            %Put on left
            leafs(node_number2).left = new_node;
         else
            leafs(node_number2).right = new_node;
         end
         
      case 2,
         %Crossover
         %First, decide where to do the crossover
         node_number1 = walk_along(solutions(s1), leafs);
         node_number2 = walk_along(solutions(s2), leafs);
         
         %Change pointers by finding which leaf points to each of the leafs to change
         l1 = []; l2 = []; r1 = []; r2 = [];
         for m = 1:length(leafs),
            if ~isempty(leafs(m).left),
               if leafs(m).left == node_number1,
                  l1 = m;
               end
               if leafs(m).left == node_number2,
                  l2 = m;
               end
            end
            if ~isempty(leafs(m).right),
               if leafs(m).right == node_number1,
	               r1 = m;
   	         end
      	      if leafs(m).right == node_number2,
         	      r2 = m;
               end
            end            
         end
         
         if ~isempty([l1 r1]) & ~isempty([l2 r2]),
            if isempty(l1),
               %It is on the right of 1
               if isempty(l2),
                  leafs(r2).right = node_number1;
                  leafs(r1).right = node_number2;
               else
                  leafs(l2).left  = node_number1;
                  leafs(r1).right = node_number2;
               end
            else
               if isempty(l2),
                  leafs(r2).right = node_number1;
                  leafs(l1).right = node_number2;
               else
                  leafs(l2).left  = node_number1;
                  leafs(l1).right = node_number2;
               end
            end
         end
         
      case 3,
         %Mutation
         %Mutate s1 or s2? 
         if rand(1) > .5,
            temp = s1;
            s1   = s2;
            s2   = temp;
         end
         
         %First, decide where to do the mutation
         node_number1 = walk_along(solutions(s1), leafs);
         
         %Choose the new operator
         o1	= 1 + floor(rand(1)*size(operators,1));
         
         %Mutate!
         leafs(node_number1).operator = operators(o1,:);
         
         if strcmp(operators(o1,1),'X'),
            leafs(node_number1).left  = [];
            leafs(node_number1).right = [];
         end
         
      end
      
   end
   
   
end


%Classify test patterns
disp(['Will use a solution with ' num2str(best*100) '% correct for the test patterns'])

dfun    = best_fun;

for k = 1:Ni,
    dfun = strrep(dfun, ['X' num2str(k)], ['test_patterns(' num2str(k) ',:)']);
end

test_targets = eval(dfun);

if (length(Uc) == 2)
    test_targets = test_targets > 0;
end

%END

function [root, leafs] = copy(start_point, leafs)
%Make a copy of a branch, from the start point until its end. Return a pointer to the copy
N	= length(leafs) + 1;
leafs(N).operator = leafs(start_point).operator;

if ~isempty(leafs(start_point).left),
   [root, leafs] = copy(leafs(start_point).left, leafs);
   leafs(N).left = root;
end
if ~isempty(leafs(start_point).right),
   [root, leafs]  = copy(leafs(start_point).right, leafs);
   leafs(N).right = root;
end

root = N;

%END


function node = walk_along(start_point, leafs);
%Find where to mutate a solution along a tree. This is recursive
r = 1 + floor(rand(1)*4);

switch r,
case {1,2},
   %Go no further
   node = start_point;
case 3,
   %If possible, go right
   if ~isempty(leafs(start_point).right),
      node = walk_along(leafs(start_point).right, leafs);
   else
      node = start_point;
   end
case 4,
   %If possible, go left
   if ~isempty(leafs(start_point).left),
      node = walk_along(leafs(start_point).left, leafs);
   else
      node = start_point;
   end
end

%END


function s = collect_tree(pointer, leafs)
%Find the tree function recursively
if ~isempty(leafs(pointer).left),
   s = ['(' collect_tree(leafs(pointer).left, leafs) ')'];
else
	s = [];   
end

s = [s leafs(pointer).operator];

if ~isempty(leafs(pointer).right),
   s = [s '(' collect_tree(leafs(pointer).right, leafs) ')'];
end

%END

function [leafs, N] = foliage(leafs, operators, remaining_depth)
%Recursively add leafs in order to build the initial tree
N			= length(leafs)+1;

if (remaining_depth > 0),
	r		= floor(rand(1)*size(operators,1)) + 1;
   leafs(N).operator = operators(r,:);
   if ~strcmp(leafs(N).operator(1), 'X'),
      %The operator is numeric, so more leafs are needed
      [leafs, added] = foliage(leafs, operators, remaining_depth-1);	      
      leafs(N).left  = added;
      [leafs, added] = foliage(leafs, operators, remaining_depth-1);	      
      leafs(N).right = added;
   else
      leafs(N).left     = [];
	   leafs(N).right    = [];
   end
else
   in = find(operators(:,1) == 'X');
	r	= floor(rand(1)*length(in)) + 1;
   leafs(N).operator = ['X' num2str(r)];
   leafs(N).left     = [];
   leafs(N).right    = [];
end
