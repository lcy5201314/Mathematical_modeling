function test_targets = Perceptron_Voted(train_patterns, train_targets, test_patterns, params)

% Classify using the Voted Perceptron algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%	Params          - [NumberOfPerceptrons, Kernel method (Linear, Polynomial, Gaussian), Kernel params]
%                     The kernel parameters are:
%			          Linear - none, Polinomial - power, Gaussian - sigma
%
% Outputs
%	test_targets	- Predicted targets
%
% NOTE: Works for only two classes
% Coded by: Igor Makienko and Victor Yosef

[NumberOfPerceptrons, method, alg_param] = process_params(params);

[c, n]		   = size(train_patterns);
train_patterns = [train_patterns ; ones(1,n)];
train_one      = find(train_targets == 1);
train_zero     = find(train_targets == 0);

%Preprocessing
processed_patterns = train_patterns;
processed_patterns(:,train_zero) = -processed_patterns(:,train_zero);

%Initial weights for Linear case:
w_percept  = rand(c+1,NumberOfPerceptrons);
%Initial alphas for kernel method:
alpha = rand(n,NumberOfPerceptrons);

%Initial permutation matrix for kernel case;

switch method
case 'Polynomial'
    perm = polyn(processed_patterns', processed_patterns',alg_param);
case 'Gaussian'
    perm = gaus(processed_patterns',processed_patterns',alg_param);
end
%Train targets for kernels' case [-1 1] :
t = 2 * train_targets - 1;
%Step for kernel case :
etta  = 1;
%Initial success vector:
w_sucesses = ones(NumberOfPerceptrons,1);

correct_classified = 0;
iter			   = 0;
max_iter		   = 500;

while (iter < max_iter)
    iter 		= iter + 1;
    indice 	= 1 + floor(rand(1)*n);
    switch method
    case 'Linear',
        InnerProduct = w_percept' * processed_patterns(:,indice);
        NegInnerProduct = (InnerProduct<=0);
        PosInnerProduct = (InnerProduct>0);
        w_sucesses = ones(size(w_sucesses)) + w_sucesses.*PosInnerProduct;
        w_percept(:,find(NegInnerProduct)) = w_percept(:,find(NegInnerProduct))...
            + processed_patterns(:,indice) * ones(1,sum(NegInnerProduct));    
        
    case {'Polynomial','Gaussian'}
        InnerProduct = perm(indice,:) * ((alpha'.*(ones(size(alpha,2),1)*t)))' ;     
        NegInnerProduct = (InnerProduct<=0)';
        PosInnerProduct = (InnerProduct>0)';
        w_sucesses = ones(size(w_sucesses)) + w_sucesses.*PosInnerProduct;
        alpha(indice,find(NegInnerProduct)) = alpha(indice,find(NegInnerProduct))...
            + etta * ones(1,sum(NegInnerProduct));
    otherwise
        error('Method unknown');
    end
end

if (iter == max_iter),
    disp(['Maximum iteration (' num2str(max_iter) ') reached'])
end

%Classify test patterns
N               = size(test_patterns, 2);
test_targets    = zeros(1, N);
test_patterns(c+1, :) = ones(1, N);

switch method
case 'Linear',
    for i = 1:NumberOfPerceptrons
        test_targets = test_targets + w_sucesses(i) * (2 * (w_percept(:,i)'*test_patterns > 0) - 1);
    end
case 'Polynomial',
    temp = [x(:),y(:),ones(size(y(:)))];
    perm  = polyn(temp,processed_patterns',alg_param);
    for i = 1:NumberOfPerceptrons,
        temp = 2 * (sum(((ones(N,1)*(alpha(:,i)'.* t)) .* (perm))') > 0) - 1;
        test_targets = test_targets + w_sucesses(i) * temp; 
    end
    
case 'Gaussian',
    temp = [x(:),y(:),ones(size(y(:)))];
    perm  = gaus(temp,processed_patterns',alg_param);
    for i = 1:NumberOfPerceptrons,
        temp = 2 * (sum(((ones(N,1)*(alpha(:,i)'.* t)) .* (perm))') > 0) - 1;
        test_targets = test_targets + w_sucesses(i) * temp; 
    end
    
end   

test_targets = test_targets > 0;

disp(['Iterated ' num2str(iter) ' times.'])

function out =  polyn(x,y,p);
%Routine function for polynomial kernel
%Input: 
%x - (number of vectors)x(dim+1) matrix  
%y - (number of vectors)x(dim+1) matrix
%p  - order of polynom
out = (ones(size(x,1),size(y,1)) + x * y').^p;

function out = gaus(x,y,sigma);
%Routine function for gaussian kernel
%Input: 
%x - (number of vectors)x(dim+1) matrix  
%y - (number of vectors)x(dim+1) matrix
%sigma  - std of gaussian kernel
x = x';y =y';c = [];
for i = 1:size(x,1),
    c(:,:,i) = (ones(size(x,2),1) * y(i,:) - x(i,:)' * ones(1,size(y,2))).^2;
end
out = exp( - squeeze( sum(permute(c,[3,1,2]))) ./ (2 * sigma) ^2);



