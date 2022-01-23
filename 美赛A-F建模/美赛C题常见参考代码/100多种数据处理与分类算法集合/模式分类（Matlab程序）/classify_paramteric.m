function targets = classify_paramteric(param_struct, patterns)

%Function for classifying patterns based on parametric distributions.
%Inputs are the paramters of the distributions and the patterns.
%Output is the classification

[Dim,N] = size(patterns);
Nc      = size(param_struct,2);
V		= zeros(Nc, N);

for j = 1:Nc,
    p = param_struct(j).p;
    for i = 1:length(param_struct(j).w),
        if (length(param_struct(j).w) > 1)
            sigma = squeeze(param_struct(j).sigma(i,:,:));
            mu    = param_struct(j).mu(i,:)';
            type  = param_struct(j).type(i);
            w     = param_struct(j).w(i);
        else
            sigma = param_struct(j).sigma;
            mu    = param_struct(j).mu';
            type  = param_struct(j).type;
            w     = param_struct(j).w;
        end
        
        if strcmp(type, 'Gaussian') 
            if length(size(sigma))>2,
                sigma = squeeze(sigma);
            end
            if (Dim == 2)
                   invsigma = inv(sigma);
                   V(j,:) = V(j,:) + w ./ (2 * pi * sqrt(abs(det(sigma)))) .* ...
                            exp(-0.5*(invsigma(1,1).*(patterns(1,:)-mu(1)).^2 + ...
                            2*invsigma(2,1).*(patterns(1,:)-mu(1)).*(patterns(2,:)-mu(2))+...
                            invsigma(2,2).*(patterns(2,:)-mu(2)).^2));
            else
                for k = 1:N,
                    V(j,k) = V(j,k) + p*w ./ ((2*pi)^(Dim/2) * sqrt(abs(det(sigma)))) .* ...
                             exp(-0.5*(patterns(:,k)-mu)'*inv(sigma)*(patterns(:,k)-mu));
                end
            end
        else
            area = prod(sigma(1,:));
            
            if (area > 0)
                min_dist = (mu-sigma(1,:)')*ones(1,N);
                max_dist = (mu+sigma(1,:)')*ones(1,N);
                
                in_range = all((patterns >= min_dist) .* (patterns <= max_dist));
                
                V(j,find(in_range)) = V(j,find(in_range)) + 1/area;
                %for k = 1:N,
                %    if (all(patterns(:,k)>mu-sigma(1,:)') & all(patterns(:,k)<mu+sigma(1,:)'))
                %        V(j,k) = V(j,k) + 1/area;
                %    end
                %end
            end
        end
  end
end

    
[m, targets]  = max(V);

%For two classes, return targets in {0,1}
if (Nc == 2)
    targets = targets - 1;
end
