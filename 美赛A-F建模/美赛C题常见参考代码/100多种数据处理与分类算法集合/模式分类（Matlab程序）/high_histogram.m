function [H, Bins, region] = high_histogram(data, Nbins, region)

%Find the histogram for high dimensional data
%
% Input:
%		Data    - The data matrix. 
%				 Should be a DxN matrix, where D is the dimension of the data, and N number of examples
%		Nbins   - Number of bins (Optional, can be a vector)
%       region  - A vector containing the min and max values in each direction
%
% Output:
%		H	  - Histogram
%       Bins  - The bin locations for the data

%Find how many bins to use
switch nargin
case 1,
   Nbins = max(3,floor(size(data,2).^(1/3)));
   region = zeros(1,2*size(data,1));
   for i = 1:size(data,1),
       region(i*2-1) = min(data(i,:));
       region(i*2)   = max(data(i,:));
   end
case 2,
   region = zeros(1,2*size(data,1));
   for i = 1:size(data,1),
       region(i*2-1) = min(data(i,:));
       region(i*2)   = max(data(i,:));
   end
case 3,
   if isempty(Nbins),
       Nbins = max(3,floor(size(data,2).^(1/3)));
   end
end

%Find the appropriate bin for each data point
Bins	= zeros(size(data));
for i = 1:size(data,1),
   one_data		= data(i,:);
   spacing		= linspace(region(i*2-1)-2*eps,region(i*2)+2*eps,Nbins(i)+1);
   minSpacing	= spacing(1:end-1);
   maxSpacing	= spacing(2:end);
   
   overMin		= zeros(Nbins(i),length(one_data));
   underMax		= zeros(Nbins(i),length(one_data));
   
   for j = 1:Nbins(i),
      overMin(j,:) 	= one_data > minSpacing(j);
      underMax(j,:)	= one_data < maxSpacing(j);      
   end
   
   [m, Bins(i,:)] = max((overMin.*underMax));
end

%Now that the data is in bins, find where they belong.
%I used three options for the dimensions, for better efficiency
switch size(data,1),
case 1,
    %1-D data
    H = zeros(Nbins,1);
    for i = 1:Nbins,
        H(i) = length(find(Bins == i));
    end
case 2,
    %2-D data
    H = zeros(Nbins);
    for i = 1:Nbins(1),
        for j = 1:Nbins(2),
            H(i,j) = length(find((Bins(1,:) == i) & (Bins(2,:) == j)));
        end
    end
otherwise
    H = zeros(Nbins*ones(1,size(data,1)));
    for i = 1:size(data,2),
        line		 = num2str(Bins(:,i))';   
        line(size(line,1)+1,:) = ones(1,size(line,2))*',';
        Fline		 = line(:)';
        eval(['H(' Fline(1:end-1) ') = H(' Fline(1:end-1) ') + 1;']);
    end
end   