function [patterns, targets, label] = SOHC(train_patterns, train_targets, Nmu, plot_on)

%Reduce the number of data points using the stepwise optimal hierarchical clustering algorithm
%Inputs:
%	train_patterns	- Input patterns
%	train_targets	- Input targets
%	Nmu 			- Number of output data points
%   plot_on         - Plot stages of the algorithm
%
%Outputs
%	patterns			- New patterns
%	targets			- New targets
%	label				- The labels given for each of the original patterns

if (nargin < 4),
    plot_on = 0;
end

[D,L]	= size(train_patterns);
c_hat   = L;

if (Nmu == 1),
   mu		= mean(train_patterns')';
   label	= ones(1,c_hat);
else  
    m       = train_patterns;
    n       = ones(1,c_hat);
    
    while (c_hat ~= Nmu),
        c_hat = c_hat - 1;
        
        %Find two clusters j and i whose merger changes the criterion the least
        %The criterion is as in DHS, Chapter 10, eq. 83
        de  = zeros(c_hat+1);
        for i = 1:c_hat+1;
            de(i,:) = sqrt(sum((m(:,i)*ones(1,c_hat+1) - m).^2)).*sqrt(n(i)*n./(n+n(i)));
            de(i,i) = 1e30;
        end
        
        [i,j] = find(de == min(min(de)));
        i     = i(1);
        j     = j(1);
        
        %Merge clusters i and j (into cluster i in this realization)
        m(:,i)  = (n(i)*m(:,i) + n(j)*m(:,j))/(n(i) + n(j));
        n(i)    = n(i) + n(j);
        if (j < c_hat + 1),
            m       = m(:,[1:j-1,j+1:end]);
            n       = n(:,[1:j-1,j+1:end]);
        else
            m       = m(:,1:end-1);
            n       = n(1:end-1);
        end
        
        if (plot_on > 0),
            plot_process(m, plot_on)
        else
            disp(['Reduced to ' num2str(c_hat) ' clusters so far'])
        end

    end
end

%Find labels for the examples
dist    = zeros(Nmu, L);
for i = 1:Nmu,
    dist(i,:) = sqrt(sum((m(:,i)*ones(1,L) - train_patterns).^2));
end
[temp, label]  = min(dist);

%Label the points
Uc        = unique(train_targets);
targets   = zeros(1,Nmu);
for i = 1:Nmu,
    N = hist(train_targets(:,find(label == i)), Uc);
    [temp, max_l] = max(N);
    targets(i) = Uc(max_l);
end

patterns = m;