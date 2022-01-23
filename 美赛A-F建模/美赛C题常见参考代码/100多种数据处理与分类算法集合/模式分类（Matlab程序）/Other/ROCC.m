function [false_alarm, hit] = ROCC(patterns, targets)

%Generate a receiver operating characteristic curve (ROCC) for 1-D data
%Inputs:
%	patterns	- The data from which to estimate
%	targets     - The class for each of the patterns
%Outputs:
%	false_alarm, hit - The x and y axes for the ROCC. 
%	If the function is called without output, the ROCC is plotted

Nbins       = max(3,floor(length(patterns).^(1/3)));
indice0		= find(targets ~= 1);
indice1		= find(targets == 1);
range		= [min(patterns), max(patterns)];

p0          = high_histogram(patterns(:,indice0),Nbins,range);
p1          = high_histogram(patterns(:,indice1),Nbins,range);

p0			= p0 ./ sum(p0);
p1			= p1 ./ sum(p1);

false_alarm = 1-cumsum(p0);
hit			= 1-cumsum(p1);

if (nargout == 0),
   figure
   plot(false_alarm, hit)
   xlabel('False alarm')
   ylabel('Hit rate')
end
