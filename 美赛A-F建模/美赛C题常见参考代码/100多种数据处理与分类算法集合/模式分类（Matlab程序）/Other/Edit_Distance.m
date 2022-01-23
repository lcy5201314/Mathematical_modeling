function C = Edit_Distance(x, y)

% Edit distance
%
% Inputs:
%	x, y				- Text vectors for the distance matrix
%
% Output:
%	C					- Distance matrix

n	= length(y);
m	= length(x);
C	= zeros(m+1,n+1);

%NOTE: Because of the matlab indexing, which starts from 1, the algorithm is slightly
%different from the pseudocode in DHS

for i = 0:m,
   C(i+1,1) = i;
end
for j = 0:n,
   C(1,j+1) = j;
end

i = 0;
j = 0;

for i = 1:m,
   for j = 1:n,
      C(i+1,j+1) = min([C(i, j+1)+1, C(i+1,j)+1, C(i,j)+1-(x(i)==y(j))]) ;
   end
end
