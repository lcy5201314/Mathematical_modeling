function test_targets = Store_Grabbag(train_patterns, train_targets, test_patterns, Knn)

% Classify using the store-grabbag algorithm (an improvement on the nearest neighbor)
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Knn		        - Number of nearest neighbors
%
% Outputs
%	test_targets	- Predicted targets

L		= length(train_patterns);

%Placing first sample in STORE
Store_patterns(:,1) = train_patterns(:,1);
Store_targets       = train_targets(1);
Grabbag_targets     = [];
Grabbag_patterns    = [];

for i = 2:L,
   target = Knn_Rule(train_patterns(:,i), Store_patterns, Store_targets, Knn);
   if target == train_targets(i)
      Grabbag_patterns = [Grabbag_patterns , train_patterns(:,i)];  
      Grabbag_targets = [Grabbag_targets train_targets(i)];
   else
      Store_patterns = [Store_patterns, train_patterns(:,i)];
      Store_targets  = [Store_targets train_targets(i)];
   end 
end      

New_Grabbag_patterns = Grabbag_patterns;

while (Grabbag_patterns ~= New_Grabbag_patterns)
   Grabbag_patterns = New_Grabbag_patterns;
   New_Grabbag_targets = [];
   for i = 1:length(Grabbag_patterns),
      target = Knn_Rule(Grabbag_patterns(:,i), Store_patterns, Store_targets);
   	if target == train_targets(i)
      	New_Grabbag_patterns = [New_Grabbag_patterns, train_patterns(:,i)];  
      	New_Grabbag_targets  = [New_Grabbag_targets train_targets(i)];
   	else
      	Store_patterns = [Store_patterns, train_patterns(:,i)];
      	Store_targets  = [Store_targets , train_targets(i)];
      end
   end
end
    
      
disp(['Calling Nearest Neighbor algorithm']);
test_targets = Nearest_Neighbor(Store_patterns, Store_targets, test_patterns, Knn);

%END

function target = Knn_Rule(Sample, Store_patterns, Store_targets, Knn)
%Classify a sample using the NN rule

for i = 1:length(Store_targets),
   %Find the k nearest neighbours
   dist(i) = sqrt((Sample(1)-Store_patterns(1,i)).^2+(Sample(2)-Store_patterns(2,i)).^2);  
end
[sorted_dist, indices] = sort(dist);

if length(Store_targets) <= Knn
   k_nearest = Store_targets;
else
   k_nearest = Store_targets(indices(1:Knn));
end

target = (sum(k_nearest) > Knn/2);

