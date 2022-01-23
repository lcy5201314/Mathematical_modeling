function [test_targets, param_struct] = EM(train_patterns, train_targets, test_patterns, Ngaussians)

% Classify using the expectation-maximization algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%   Ngaussians      - Number for Gaussians for each class (vector)
%
% Outputs
%	test_targets	- Predicted targets
%   param_struct    - A parameter structure containing the parameters of the Gaussians found

classes             = unique(train_targets); %Number of classes in targets
Nclasses            = length(classes);
Nalpha				= Ngaussians;						 %Number of Gaussians in each class
Dim                 = size(train_patterns,1);

max_iter   			= 100;
max_try             = 5;
Pw					= zeros(Nclasses,max(Ngaussians));
sigma				= zeros(Nclasses,max(Ngaussians),size(train_patterns,1),size(train_patterns,1));
mu					= zeros(Nclasses,max(Ngaussians),size(train_patterns,1));

%The initial guess is based on k-means preprocessing. If it does not converge after
%max_iter iterations, a random guess is used.
disp('Using k-means for initial guess')
for i = 1:Nclasses,
    in  			= find(train_targets==classes(i));
    [initial_mu, targets, labels]	= k_means(train_patterns(:,in),train_targets(:,in),Ngaussians(i));
    for j = 1:Ngaussians(i),
        gauss_labels    = find(labels==j);
        Pw(i,j)         = length(gauss_labels) / length(labels);
        sigma(i,j,:,:)  = diag(std(train_patterns(:,in(gauss_labels))'));
    end
    mu(i,1:Ngaussians(i),:) = initial_mu';
end

%Do the EM: Estimate mean and covariance for each class 
for c = 1:Nclasses,
    train	   = find(train_targets == classes(c));
    
    if (Ngaussians(c) == 1),
        %If there is only one Gaussian, there is no need to do a whole EM procedure
        sigma(c,1,:,:)  = sqrtm(cov(train_patterns(:,train)',1));
        mu(c,1,:)       = mean(train_patterns(:,train)');
    else
        
        sigma_i         = squeeze(sigma(c,:,:,:));
        old_sigma       = zeros(size(sigma_i)); 		%Used for the stopping criterion
        iter			= 0;									%Iteration counter
        n			 	= length(train);					%Number of training points
        qi			    = zeros(Nalpha(c),n);	   	%This will hold qi's
        P				= zeros(1,Nalpha(c));
        Ntry            = 0;
        
        while ((sum(sum(sum(abs(sigma_i-old_sigma)))) > 1e-4) & (Ntry < max_try))
            old_sigma = sigma_i;
            
            %E step: Compute Q(theta; theta_i)
            for t = 1:n,
                data  = train_patterns(:,train(t));
                for k = 1:Nalpha(c),
                    P(k) = Pw(c,k) * p_single(data, squeeze(mu(c,k,:)), squeeze(sigma_i(k,:,:)));
                end          
                
                for i = 1:Nalpha(c),
                    qi(i,t) = P(i) / sum(P);
                end
            end
            
            %M step: theta_i+1 <- argmax(Q(theta; theta_i))
            %In the implementation given here, the goal is to find the distribution of the Gaussians using
            %maximum likelihod estimation, as shown in section 10.4.2 of DHS
            
            %Calculating mu's
            for i = 1:Nalpha(c),
                mu(c,i,:) = sum((train_patterns(:,train).*(ones(Dim,1)*qi(i,:)))')/sum(qi(i,:)');
            end
            
            %Calculating sigma's
            %A bit different from the handouts, but much more efficient
            for i = 1:Nalpha(c),
                data_vec = train_patterns(:,train);
                data_vec = data_vec - squeeze(mu(c,i,:)) * ones(1,n);
                data_vec = data_vec .* (ones(Dim,1) * sqrt(qi(i,:)));
                sigma_i(i,:,:) = sqrt(abs(cov(data_vec',1)*n/sum(qi(i,:)')));
            end
            
            %Calculating alpha's
            Pw(c,1:Ngaussians(c)) = 1/n*sum(qi');
            
            iter = iter + 1;
            disp(['Iteration: ' num2str(iter)])
            
            if (iter > max_iter),
                theta = randn(size(sigma_i));
                iter  = 0;
                Ntry  = Ntry + 1;
                
                if (Ntry > max_try)
                    disp(['Could not converge after ' num2str(Ntry-2) ' redraws. Quitting']);
                else
                    disp('Redrawing weights.')
                end
            end
            
        end
        
        sigma(c,:,:,:) = sigma_i;
    end
end

%Classify test patterns
for c = 1:Nclasses,
    param_struct(c).p       = length(find(train_targets == classes(c)))/length(train_targets);
    param_struct(c).mu      = squeeze(mu(c,1:Ngaussians(c),:));
    param_struct(c).sigma   = squeeze(sigma(c,1:Ngaussians(c),:,:));
    param_struct(c).w       = Pw(c,1:Ngaussians(c));
    for j = 1:Ngaussians(c)
        param_struct(c).type(j,:) = cellstr('Gaussian');
    end
    if (Ngaussians(c) == 1)
        param_struct(c).mu = param_struct(c).mu';
    end
end
test_targets = classify_paramteric(param_struct, test_patterns);

%END EM

function p = p_single(x, mu, sigma)

%Return the probability on a Gaussian probability function. Used by EM

p = 1/(2*pi*abs(det(sigma)))^(length(mu)/2)*exp(-0.5*(x-mu)'*inv(sigma)*(x-mu));

