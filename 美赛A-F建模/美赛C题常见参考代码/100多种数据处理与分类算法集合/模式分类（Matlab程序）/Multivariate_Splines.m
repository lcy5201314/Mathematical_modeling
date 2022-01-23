function test_targets = Multivariate_Splines(train_patterns, train_targets, test_patterns, params)

% Classify using multivariate adaptive regression splines
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	params          - [Spline degree, Number of knots per spline]
%
% Outputs
%	test_targets	- Predicted targets

[Ni, M]		    = size(train_patterns);
options         = optimset('Display','off');

%Get parameters
[q, Nk]			= process_params(params);

order           = zeros(1,Ni); %The new order of the dimensions
residual        = train_targets;
knots           = zeros(Nk, Ni);

%Define what a spline is
Pspline         = inline('sum(((x-t*ones(1,L)).^q).*(((x-t*ones(1,L)).^q)>0))', 't', 'x', 'q', 'L');
minPspline      = inline('sum((targets - sum(((x-t*ones(1,L)).^q).*(((x-t*ones(1,L)).^q)>0))).^2)/L', 't', 'x', 'q', 'L', 'targets');

for i = 1:Ni,
    %Find which remaining dimension is fit by a spline:
    
    %Findd the remaining dimensions
    remaining = find(~ismember(1:Ni, order));
    
    %Fit a spline to each of the dimensions
    temp_knots  = zeros(Nk, Ni-i+1);
    errors      = zeros(1, Ni-i+1);
    for j = 1:Ni-i+1,
        temp_knots(:,j) = fminunc(minPspline, randn(Nk, 1), options, ones(Nk,1)*train_patterns(remaining(j),:), q, M, residual);
        errors(j)       = feval(minPspline, temp_knots(:,j), ones(Nk,1)*train_patterns(remaining(j),:), q, M, residual);
    end
    
    %Find the best dimension to regress on
    [best, best_dim] = min(errors);
    order(i)         = remaining(best_dim);
    knots(:,i)       = temp_knots(:,best_dim);
    
    %Compute residual 
    predict          = feval(Pspline, temp_knots(:,j), ones(Nk,1)*train_patterns(remaining(j),:), q, M);
    residual         = residual ./ predict;
end

%Compute weights via pseudo-inverse:
%Compute the prediction for each dimension
prediction = zeros(Ni, M);
for i = 1:Ni,
    prediction(i,:) = feval(Pspline, knots(:,i), ones(Nk,1)*train_patterns(order(i),:), q, M);
    if i > 1,
        prediction(i,:) = prod(prediction(1:i,:));
    end
end

%Compute the weights
W   = pinv(prediction*prediction')*prediction*train_targets';

%Compute the test targets
N           = size(test_patterns,2);
prediction  = zeros(Ni, N);
for i = 1:Ni,
    prediction(i,:) = feval(Pspline, knots(:,i), ones(Nk,1)*test_patterns(order(i),:), q, N);
    if i > 1,
        prediction(i,:) = prod(prediction(1:i,:));
    end
end

test_targets = W'*prediction;

if (length(unique(train_targets)) == 2)
    test_targets = test_targets > 0.5;
end
