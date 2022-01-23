function [test_targets, P, theta, phi] = LocBoost(train_patterns, train_targets, test_patterns, params)

% Classify using the local boosting algorithm
% Inputs:
% 	train_patterns	- Train patterns
%	train_targets	- Train targets
%   test_patterns   - Test  patterns
%   params          - A vector containing the algorithm paramters:
%                     [Number of boosting iterations, number of EM iterations, Number of optimization iterations, Weak learner, Weak learner parameters]
%                     IMPORTANT: The weak learner must return a hyperplane parameter vector, as in LS
%
% Outputs
%	test_targets	- Predicted targets
%   P               - The probability function (NOT the probability for the train_targets!)
%	theta		    - Sub-classifier parameters
%	phi		        - Sub-classifier weights

th = 0.6;
test_percentage = 0.1;                  %Percentage of points to be used as a test set
[Dims, Nf]      = size(train_patterns);
Nt              = 0;
train_targets	= (train_targets > .5)*2-1;	%Set targets to {-1,1}
[Niterations, Nem, Noptim, Wtype, Wparams] = process_params(params);
Niterations		= Niterations + 1;
dist            = [];
errors          = ones(1, Niterations);

%if ((Niterations < 1) | (Nem < 1) | (Noptim < 1)),
%   error('Iteration paramters must be positive!');
%end

options1		= optimset('Display', 'off', 'LevenbergMarquardt','on');
options2		= optimset('Display', 'off', 'MaxIter', Noptim,'LevenbergMarquardt','on');

%Find first iteration parameters
theta     	= zeros(1, Dims+1);
phi			= zeros(1, Dims+Dims^2);
h			= ones(1,Nf);
counter     = 1;

%Initial value is the largest connected component again all the others
[D, tmp_theta]                 = feval(Wtype, train_patterns, train_targets, train_patterns(:,1:2), Wparams);
theta(1, 1:size(tmp_theta, 2)) = tmp_theta;

P           = 0.5*ones(1,Nf); 
%P           = LocBoostFunctions(theta(1,:), 'class_kernel', [train_patterns; ones(1,Nf)], train_targets); 
errors(1)   = mean(P<.5);
max_log_likelihood = sum(log(P));
stop        = 0;

%Find the local classifiers
for t = 2:Niterations,
    if (stop == 1)
        break;
    end
    
    %Do inital guesses
    recompute_initial_values = 0;
    [components, dist]	= compute_initial_value(train_patterns, (P<th), dist);
    if (all(components == 0))
        phi     = phi(1:counter,:);
        theta   = theta(1:counter,:);
        break
    end
        
    Uc                  = unique(components);
    Nc                  = hist(components, Uc);
    in                  = find(Uc>0);
    
    Uc                  = Uc(in);
    Nc                  = Nc(in);
    
    [Nc, si] = sort(Nc);
    Uc       = Uc(si);
    
    Uc       = fliplr(Uc);
    Nc       = fliplr(Nc);
    
    length_one  = sum(Nc == 1);
%     if (all(Nc == 1))
%         all_one = 1;
%     else
%         all_one = 0;
%     end
    
    if isempty(components) %| all_one
        phi     = phi(1:counter,:);
        theta   = theta(1:counter,:);
        break
    end
    
%     if (all_one & (errors(counter) > errors(counter-1)))
%        phi         = phi(1:counter-1,:);
%        theta       = theta(1:counter-1,:);
%        break
%     end
    
%     if (all_one == 1)
%         Uc = 1;
%     end
%     
    for i = 1:length(Uc),
        if (length_one ~= sum(Nc == 1))
            indices = find(components == Uc(i));
        else
            rnd     = randperm(length(Uc));
            indices = find(components == Uc(rnd(1)));
           
            
%             %Find the two closest samples
%             new_dist = dist(indices,indices);
%             new_dist = new_dist + eye(length(in)).*1e10;
%             [r,c]    = find(min(min(new_dist)) == new_dist);
%             indices  = indices([r(1), c(1)]);
%             %[p, t, labels] = k_means(train_patterns(:,indices), train_targets(indices), floor(length(in)/10), 0);
%             %Uc             = unique(labels);
%             %Nc             = hist(labels, Uc);
%             %[m, best]      = max(Nc);
%             %indices        = indices(find(labels == best));
        end
        
        %plot_process(train_patterns(:,indices),1)
        
        if ((length_one ~= sum(Nc)) & (length(indices) == 1))
            continue;
        end
        
        counter = counter + 1;
        
        No_EM = 0;
        
        if (length(indices) > 100)
            disp('Using k-means')
            try
                [p, t, labels] = k_means(train_patterns(:,indices), train_targets(indices), 2, 0);
                Ucl            = unique(labels);
                Ncl            = hist(labels, Ucl);
                [m, best]      = max(Ncl);
                indices        = indices(find(labels == Ucl(best)));            
                recompute_initial_values = 1;
            catch
            end
        end
        
        if (length(indices) > 1)
            means          = mean(train_patterns(:,indices)');
            stds           = var(train_patterns(:,indices)');
            covar          = cov(train_patterns(:,indices)',1)';
            
            if (all(abs(stds)<1e-5))
                full_stds      = eye(Dims)*1e3;
                No_EM = 1;
            else
                if (cond(covar) < 1e10)
                   full_stds      = inv(covar); 
               else
                   if any(abs(stds)<1e-5)
                        zs         = find(abs(stds)<1e-5);
                        temp       = stds;
                        temp(zs)   = 1;
                        temp       = 1./temp.^2;
                        temp(zs)   = 1e5;
                        full_stds  = diag(temp);
                    else
                        full_stds  = inv(diag(stds.^2));
                    end
                end
            end
        else
            means          = train_patterns(:,indices)';
            full_stds      = eye(Dims)*1e3;
            No_EM = 1;
        end

        phi(counter,1:Dims)     = means;
        phi(counter,Dims+1:end) = full_stds(:)';

        if (No_EM == 1) & is_another_example_from_other_class_close_by(train_patterns, train_targets, indices)
            length_one = length_one - 1;
            counter = counter - 1;
            phi     = phi(1:counter,:);
            theta   = theta(1:counter,:);
            continue;
        end      

        if (length(unique(train_targets(indices))) > 1)
            [D, theta(counter, 1:size(theta, 2))] = feval(Wtype, train_patterns(:,indices), train_targets(indices), train_patterns(:,1:2), Wparams);
            No_EM = 0;
        else
            theta(counter,:)   = 0;
            theta(counter,end) = unique(train_targets(indices))*10;
        end

        if (No_EM == 0)
            for i = 1:Nem,
                %Compute h(t-1)
                gamma_ker 	   = real(LocBoostFunctions(phi(counter,:), 'gamma_kernel', train_patterns, [], [], Dims));  %Gamma(x, gamma(C))
                class_ker      = LocBoostFunctions(theta(counter,:), 'class_kernel', [train_patterns; ones(1,Nf)], train_targets); 
                h_tminus1      = gamma_ker .* class_ker ./ ((1-gamma_ker).*P + gamma_ker.*class_ker);
	
                %Optimize theta(t,:) using first part of the Q function
                if (length(unique(train_targets(indices))) > 1)
                    temp_theta     = fminsearch('LocBoostFunctions', theta(counter,:), options1, 'Q1', [train_patterns; ones(1,Nf)], train_targets, h_tminus1);
                    %[D, temp_theta] = feval(Wtype, train_patterns, train_targets, train_patterns(:,1:2), h_tminus1);
                else
                    temp_theta      = zeros(1,Dims+1);
                    temp_theta(end) = unique(train_targets(indices))*10;
                end         
                %[d, temp_theta(1,1:size(theta,2))] = feval('LS', train_patterns, train_targets, train_patterns(:,1:2), h_tminus1);
                
                %Optimize gamma(t,:) using second part of the Q function
                %temp_phi       = fmincon('LocBoostFunctions', phi(t,:), [], [], [], [], lb, [], [], options, 'Q2',  train_patterns, train_targets, h_tminus1, Dims);
                temp_phi       = fminsearch('LocBoostFunctions', phi(counter,:), options2, 'Q2',  train_patterns, train_targets, h_tminus1, Dims);
                
                theta(counter,:) = temp_theta;
                phi(counter,:)   = temp_phi;
            end
        end 
        
        oldP = P;
        %Compute new P function 
        gamma_ker	   = real(LocBoostFunctions(phi(counter,:), 'gamma_kernel', train_patterns, [], [], Dims));  
        class_ker      = LocBoostFunctions(theta(counter,:), 'class_kernel', [train_patterns; ones(1,Nf)], train_targets); 
        P              = max(eps, (1-gamma_ker).*P + gamma_ker.*class_ker);
        log_likelihood = sum(log(P));
        
        errors(counter)     = mean(P<.5);
        
        %figure(2)
        %contourf(reshape(LocBoostFunctions(phi(1:counter,:), 'NewTestSet', test_patterns, ones(1, size(test_patterns,2)), [], theta(1:counter,:))>.5,100,100))
        %figure(1)
        disp(['Iteration ' num2str(counter-1) ': Train error: ' num2str(sum(P<.5)/length(P)) ', Log likelihood ' num2str(log_likelihood) ', Component length ' num2str(length(indices))])
        
        if (((length_one == sum(Nc == 1)) & (errors(counter) > errors(counter-1))) | (recompute_initial_values == 1))
            break
        end
        
        if ((log_likelihood < max_log_likelihood - 0) & (errors(counter) > errors(counter-1))),
            %Nothing more to do
            phi         = phi(1:counter-1,:);
            theta       = theta(1:counter-1,:);
            disp('Jump in Log likelihood')
            stop        = 1;
            break
        end
 
        if (log_likelihood == 0),
            %Nothing more to do
            phi         = phi(1:counter,:);
            theta       = theta(1:counter,:);
            stop        = 1;
            disp('Log likelihood is zero')
            break
        end
        
        if (counter > Niterations)
            phi         = phi(1:counter,:);
            theta       = theta(1:counter,:);
            stop        = 1;
            disp('Max iteration reached')
            break
        end
            
        if (~isreal(log_likelihood) | ~isfinite(log_likelihood)),
            %Nothing more to do
            phi         = phi(1:counter-1,:);
            theta       = theta(1:counter-1,:);
            stop        = 1;
            disp('Log likelihood is not finite or is not real')
            break
        end
        
        if (max_log_likelihood < log_likelihood)
            max_log_likelihood = log_likelihood;
        end
    end
    
end

% [m, cut] = min(errors); cut = cut(1);
% phi      = phi(1:cut,:);
% theta    = theta(1:cut,:);

%Classify test patterns
test_targets = LocBoostFunctions(phi, 'NewTestSet', test_patterns, ones(1, size(test_patterns,2)), [], theta);

test_targets = test_targets > 0.5;

%end LocBoost
%*********************************************************************

function [component, dist] = compute_initial_value(train_patterns, train_targets, dist)

%Returns the initial guess by connected components

[Dim,n] = size(train_patterns);

% Compute all distances, if it has not been done before
if (isempty(dist)),
    dist = zeros(n);
    for i = 1:n,
        dist(i,:) = sum((train_patterns(:,i)*ones(1,n) - train_patterns).^2);
    end
end

ind_plus	= find(train_targets == 1);
size_plus   = length(ind_plus);

G = zeros(n);
for i=1:size_plus   
    [o,I] = sort(dist(ind_plus(i),:));
    for j=1:n
        if (train_targets(I(j)) == 1),
            G(ind_plus(i),I(j)) = 1;
            G(I(j),ind_plus(i)) = 1;
        else
            break
        end
    end
end
G = G - (tril(G).*triu(G)); %Remove main diagonal

if ~all(diag(G)) 
    [p,p,r,r] = dmperm(G|speye(size(G)));
else
    [p,p,r,r] = dmperm(G);  
end;

% Now the i-th component of G(p,p) is r(i):r(i+1)-1.
sizes   = diff(r);        % Sizes of components, in vertices.
k       = length(sizes);      % Number of components.

% Now compute an array "blocks" that maps vertices of G to components;
% First, it will map vertices of G(p,p) to components...
component           = zeros(1,n);
component(r(1:k))   = ones(1,k);
component           = cumsum(component);

% Second, permute it so it maps vertices of A to components.
component(p) = component;

component    = component .* train_targets; %Mark all correctly assigned targets as zeros

function ret = is_another_example_from_other_class_close_by (train_patterns, train_targets, indices)

if (length(unique(train_targets(indices))) ~= 1)
    ret = 1;
    return;
end

buffer = indices;

for i = 1:size(train_patterns,2)
    if sum(ismember(buffer, i) == 0)
        dist = sum((train_patterns(:,i) - train_patterns(:,indices(1))).^2);
        if (dist < 0.1)
            buffer = [buffer, i];
        end
    end
end

if (length(unique(train_targets(buffer))) ~= 1)
    ret = 1;
else
    ret = 0;
end

        