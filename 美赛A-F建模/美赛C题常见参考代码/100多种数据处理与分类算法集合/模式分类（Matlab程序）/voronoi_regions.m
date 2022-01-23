function D = voronoi_regions(patterns, region)

% Make a Voronoi diagram from sample points
% Inputs:
%	patterns	- Input data patterns
%	targets	- Input data targets
%	region	- Decision region vector: [-x x -y y number_of_points]

N		= region(5);
x		= linspace (region(1),region(2),N);
y		= linspace (region(3),region(4),N);
D		= zeros(N);
[r,c] = size(patterns);

y_dist	= (ones(N,1) * patterns(2,:) - y'*ones(1,c)).^2;
for i = 1:N,
	if (i/50 == floor(i/50)),
      disp(['Finished ' num2str(i) ' lines out of ' num2str(N) ' lines.'])
   end
   x_dist = ones(N,1)  * (patterns(1,:)-x(i)).^2;
   dist   = abs(x_dist + y_dist);   
   [sorted_dist, indices] = min(dist');
   D(:,i) = indices(1,:)';
end