function [patterns, targets] = generate_data_set(parameters)

% Generate a new data set given it's Gaussian parameters
% Inputs:
%	parameters: A parameteric distribution structure
%	region	- Decision region vector: [-x x -y y number_of_points]
% 
% Outputs:
%	patterns - New patterns
%	targets  - New targets 

%Get the number of points from the user
N       = str2num(char(inputdlg('Enter the number of points in data set:','Generate Data set')));
Uc      = length(parameters);
Dim     = size(parameters(1).mu,2);
targets = zeros(1, N);
patterns= zeros(Dim, N);

Nv  = zeros(1, Uc);
for i = 1:Uc,
    Nv(i) = parameters(i).p;
end
Nv = round([0, cumsum(Nv)]*N);

%First, select the indices for the different classes
for j = 1:length(parameters),
    indices             = [Nv(j)+1:Nv(j+1)];
    targets(indices)    = j-1;
    
    %Now, make the points for this class
    %First, divide the indices into pieces according to w
    w       = round(length(indices)*parameters(j).w);
    cw      = [0; cumsum(w)];
    cw(end) = length(indices);
    
    for i = 1:length(w),
        %Make w(i) data points according to the distribution
        if (length(size(parameters(j).sigma))>2),
            sigma = squeeze(parameters(j).sigma(i,:,:));
        else
            sigma = parameters(j).sigma;
        end
        
        if (strcmp('Gaussian', parameters(j).type(i)))
            A		= sigma;
            dist 	= A'*randn(Dim,w(i));
            dist    = dist + parameters(j).mu(i,:)' * ones(1,w(i));
        else
            %Uniform distribution
            A       = sigma(1,:)';
            dist    = (A*ones(1,w(i))).*rand(Dim,w(i)) + (parameters(j).mu(i,:)'-A/2) * ones(1,w(i));
        end
        
        %Place them in one of the remaining places of indices0
	    patterns(:,indices(1+cw(i):cw(i+1))) = dist(:,1:cw(i+1)-cw(i));
    end
end

%Mix the data
indices  = randperm(N);
patterns = patterns(:,indices);
targets  = targets(indices);