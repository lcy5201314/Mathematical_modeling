function [new_patterns, targets] = MDS(patterns, targets, params)

% Reshape the data using the multidimensional scaling algorithm
% Inputs:
% 	patterns	- Train patterns
%	targets	- Train targets
%	params   - [Method, Number of output dimensions, Convergence rate]
%
% Outputs
%	patterns - New patterns
%	targets  - Targets
%
% The method can be one of three:
%	ee		- Emphesize large errors
%	ef		- Emphesize large products of errors and fractional errors
%	ff		- Emphesize large fractional errors


maxIter = 1000;
iter    = 0;
[Ni, M] = size(patterns);

[method, No, eta] = process_params(params);

if (No > Ni),
   error('Error: The output dimension should be smaller than the input dimension')
end

y			= PCA(patterns, targets, No);

AbsGradJ = inf;
gradJ		= zeros(size(y));

%Compute the distance matrix in the original dimension
temp     = repmat(patterns,[1 1 M]);
d		   = sqrt(squeeze(sum((temp - permute(temp, [1 3 2])).^2)));
   
while (AbsGradJ > 0.01),
   %Compute the new delta (distances in the new dimension)
	temp     = repmat(y,[1 1 M]);
	dif	   = temp - permute(temp, [1 3 2]);
	delta	   = sqrt(squeeze(sum(dif).^2));
   
   %Compute the gradient
   for k = 1:M,
      indices		= [1:k-1,k+1:M];
      switch method,
      case 'ee',
         sigma			= ones(No,1)*((d(k,indices) - delta(k,indices))./d(k, indices)) .* ...
         				  (y(:,k)*ones(1,M-1) - y(:,indices));
         gradJ(:,k)	= eta * (2/sum(sum(triu(delta,1).^2))*sum(sigma'))';
      case 'ef',
         sigma			= ones(No,1)*((d(k,indices) - delta(k,indices))./(d(k, indices).*delta(k, indices).^2)) .* ...
         				  (y(:,k)*ones(1,M-1) - y(:,indices));
         gradJ(:,k)	= eta * 2*sum(sigma')';
      case 'ff',              
         sigma			= ones(No,1)*((d(k,indices) - delta(k,indices))./(d(k, indices).*delta(k, indices).^2)) .* ...
         				  (y(:,k)*ones(1,M-1) - y(:,indices));
         gradJ(:,k)	= eta * (2/sum(sum(triu(delta,1)))*sum(sigma'))';
      otherwise
         error('Unknown method')
      end
      
   end
   
   if (AbsGradJ*2 < sum(sum(abs(gradJ)))),
      %The update got much larger
      break
   end
   
   %Update y
   y			= y - gradJ;
   
   AbsGradJ = sum(sum(abs(gradJ)));
   iter		= iter + 1;
   if (iter / 50 == floor(iter/50)),
      disp(['Iteration ' num2str(iter) ': Update is ' num2str(AbsGradJ)])
   end
   
end   

disp(['Converged after ' num2str(iter) ' iterations'])

new_patterns = y;