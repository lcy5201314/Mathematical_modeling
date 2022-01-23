function [classify, err] = classification_error(D, patterns, targets, region)

%Find a classification error for a given decision surface D and a given set of
%patterns (2xL) and targets (1xL)
%The patterns are assumed to lay in a box defined by region
%([-x x -y y ... number_of_points])


N		= region(length(region));			%Number of points on the grid
N0		= length(find(targets == 0));		%Number of targets from class 0
N1		= length(find(targets == 1));		%Number of targets from class 1
Dims	= size(patterns,1);					%Number of data dimentions
ok0	= 0;
ok1	= 0;

%Normalize the patterns according to the grid, and then check each point if it is
%classified correctly or not.
for i = 1:Dims,
   base = (i-1)*2;
	patterns(i,:) = round(1/(region(base+2)-region(base+1))*...
      ((N-1)*patterns(i,:)+region(base+2)-N*region(base+1)));
end

switch Dims,
case 2,
   for i = 1:N0+N1,
      %Get the pattern, but if for some reason it is out of range, put it on the edge of the range
      pattern1 = min(max(1,patterns(1,i)),region(5));
      pattern2 = min(max(1,patterns(2,i)),region(5));
   	if (targets(i)),
      	ok1		 = ok1 + (D(pattern2,pattern1) == 1);
	   else    
   	   ok0		 = ok0 + (D(pattern2,pattern1) == 0);
      end 
   end
case 3,
   for i = 1:N0+N1,
      pattern1 = min(max(1,patterns(1,i)),region(5));
      pattern2 = min(max(1,patterns(2,i)),region(5));
      pattern3 = min(max(1,patterns(3,i)),region(5));
   	if (targets(i)),
      	ok1		 = ok1 + (D(pattern2,pattern1,pattern3) == 1);
	   else    
   	   ok0		 = ok0 + (D(pattern2,pattern1,pattern3) == 0);
      end 
   end
otherwise
   error('The toolbox can assess the error for two or three dimensional data only.')
end


classify 		= zeros(2);
if (N0 > 0),
    classify(1,1)  = ok0/N0;
    classify(1,2)  = (N0-ok0)/N0;
end
if (N1 > 0),
    classify(2,2)  = ok1/N1;
    classify(2,1)  = (N1-ok1)/N1;
end

err = (N0*classify(1,2) + N1*classify(2,1))/(N0+N1);
