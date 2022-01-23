function test_targets = PNN(train_patterns, train_targets, test_patterns, sigma)

% Classify using a probabilistic neural network
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	sigma           - Gaussian width
%
% Outputs
%	test_targets	- Predicted targets

[Dim, Nf]       = size(train_patterns);
Dim             = Dim + 1;
train_patterns(Dim,:) = ones(1,Nf);
u_targets       = unique(train_targets);

%Build the classifier
x               = train_patterns;
W               = x ./ (ones(Dim,1)*sqrt(sum(x.^2)));  %x_jk <- x_jk / sqrt(sum(x_ji^2)), w_jk <- x_jk

%if x in w_i then a_ji <- 1
a = zeros(Nf, length(u_targets));
for i = 1:length(u_targets),
    a(find(train_targets == u_targets(i)),i) = 1;
end

%Test it and classify the test patterns
test_patterns = [test_patterns; ones(1, size(test_patterns,2))];
test_patterns = test_patterns ./ (ones(Dim,1)*sqrt(sum(test_patterns.^2)));

%net_k <- W'_t*x
net         = W' * test_patterns;

%if a_ki=1 then g_i <- g_i + exp((net-1)/sigma^2)
arguments   = zeros(length(u_targets),size(test_patterns,2));
for i = 1:length(u_targets),
    mask            = a(:,i) * ones(1,size(test_patterns,2));
    arguments(i,:)  = sum(exp((net-1)/sigma^2) .* mask);
end

%class <- argmax g(x)
[m, indices] = max(arguments);
test_targets = zeros(1,size(test_patterns,2));
for i = 1:length(u_targets),
    test_targets(find(indices == i)) = u_targets(i);
end


    