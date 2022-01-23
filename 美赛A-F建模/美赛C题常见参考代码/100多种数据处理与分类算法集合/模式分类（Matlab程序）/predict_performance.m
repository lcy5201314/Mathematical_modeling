function a = predict_performance(algorithm, algorithm_params, patterns, targets)

% Predict the final performance of an algorithm from the learning curves
% Inputs:
%   algorithm           - The algorithm to test
%   algorithm_params    - Algorithm parameters
% 	patterns            - Train patterns
%	targets	            - Train targets
%
% Outputs
%	a		        	- Final performance prediction

[Ni, M] = size(patterns);

%Define the points in the data to check
Npoints = 10;
Nfold   = 5;
points  = linspace(0.05, 0.9, Npoints);

Etest   = zeros(Nfold, Npoints);
Etrain  = zeros(Nfold, Npoints);

%Train and test the classifiers
for i = 1:Npoints,
    for j = 1:Nfold,
        in      = randperm(M);
        train_in= in(1:floor(M*points(i)));
        test_in = in(floor(0.9*M)+1:end);
        predicted_targets = feval(algorithm, patterns(:, train_in), targets(train_in), [patterns(:, train_in), patterns(:, test_in)], algorithm_params);

        etrain(j,i)      = mean(targets(train_in) ~= predicted_targets(1:length(train_in)));     
        etest(j,i)       = mean(targets(test_in) ~=  predicted_targets(1+length(train_in):end));     
    end
    disp(['Finished ' num2str(i/Npoints*100) '% for this algorithm'])
        
end

%Find the parameters of the distribution
Etest   = mean(etest);
Etrain  = mean(etrain);

%alpha is the slope of the difference curve on a log-log scale, so:
p       = polyfit(log(floor(M*points)), log(Etest-Etrain),1);
alpha   = -p(2);

A       = [ones(Npoints,1) floor(M*points').^(-alpha) floor(M*points').^(-alpha)];
b       = (Etest+Etrain)';
X       = pinv(A)*b;

a       = real(X(1)/2);