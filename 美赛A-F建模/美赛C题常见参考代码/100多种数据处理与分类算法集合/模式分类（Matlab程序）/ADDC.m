function [patterns, targets] = ADDC(train_patterns, train_targets, Nmu, plot_on)

%Reduce the number of data points using the Agglomerative clustering algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	Nmu				- Maximum number of output data points
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns		- New patterns
%	targets			- New targets

if (nargin < 4),
    plot_on = 0;
end

if (Nmu == 1),
    %If one center is needed, it is simply the average of the data
    patterns = mean(train_patterns')';
    targets  = (sum(train_targets)/length(train_targets) > 0.5);
    return
end
    
[D,L]				= size(train_patterns);
min_percentage = 0.001; %Points with count less than this will be removed
min_number		= 5; 		%Points with count less than this will also be removed

%Initialize the mu's
K		= 0; %Number of centroids
mu		= zeros(D,Nmu);
count   = zeros(1,Nmu);
Uc      = unique(train_targets);

for i = 1:L,
   data = train_patterns(:,i);
   
   if (K > 0),
      %Find closest centriod
      dist 				= sum((mu(:,1:K) - data * ones(1,K)).^2);
      [temp, min_d]  = min(dist);
      mu(:,min_d)		= mu(:,min_d) + (data - mu(:,min_d)) / (count(:,min_d) + 1);
      count(:,min_d) = count(:,min_d) + 1;
   end
   
   if (K < Nmu),
	   %Add new centroid
      K = K + 1;
      mu(:,K) = data;
   else
      %Merge redundant centroids
      closest_i1 = 0;
      closest_i2 = 0;
      dist		  = 1e100;
      for i1 = 1:K,
         for i2 = 1:K,
            if (i1 ~= i2),
               temp_dist = norm(mu(:,i1)-mu(:,i2));
               if (temp_dist < dist),
                  dist = temp_dist;
                  closest_i1 = i1;
                  closest_i2 = i2;
               end
            end
         end
      end
      if ((count(closest_i1) + count(closest_i2)) > 0),
         mu(:,closest_i1)  = (mu(:,closest_i1)*count(closest_i1) + mu(:,closest_i2)*count(closest_i2)) / ...
      						     (count(closest_i1) + count(closest_i2));
	      count(closest_i1) = count(closest_i1) + count(closest_i2);
   	   mu(:,closest_i2)  = data;
         count(closest_i2) = 0;
      end
   end
   
      %Plot the centers during the process 
      plot_process(mu, plot_on)
   
end

%Post-processing
keep 		= find(count(1:K) > max(min_percentage*L,min_number));
patterns = mu(:,keep);
Nmu		= length(keep);

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
