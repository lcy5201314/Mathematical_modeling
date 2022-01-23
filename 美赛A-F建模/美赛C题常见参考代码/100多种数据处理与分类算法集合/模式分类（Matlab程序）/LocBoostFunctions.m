function [f, Df] = LocBoostFunctions(params, type, patterns, targets, h, params2)

%Return a value for the LocBoost algorithm functions

[c,r] = size(params);
r2		= r/2;
Nf		= size(patterns,2);

switch type
case 'class_kernel'
    w   = 1./(1 + exp(-params* patterns));
	f   = (w.^((1+targets)/2) .* ((1-w).^((1-targets)/2)));
case 'Q1'
   f = LocBoostFunctions(params, 'class_kernel', patterns, targets);

   f = -sum(h.*log(eps + f));
   %f = -sum(h.*((1+targets)/2.*log(eps+(1./(1 + exp((-params* patterns))))) + ...
   %   			 (1-targets)/2.*log(eps+1-(1./(1 + exp((-params* patterns)))))));
case 'gamma_kernel'
	%f = (1./(1 + exp(params(1:r2) * patterns)));
    mu          = params(1:params2)';
    invsigma    = sqrtm(reshape(params(1+params2:end), params2, params2));

    patterns    = patterns - mu*ones(1,Nf);
    
    if (Nf < 50)
        f       = diag(patterns'*invsigma*patterns)';
    else
        f = zeros(1,Nf);
        for i = 1:Nf,
            f(i) = patterns(:,i)'*invsigma*patterns(:,i);
        end
    end
    
    f = exp(-0.5*f);
    %f = f < 1;
    %f       = exp(-0.5*sum((params(1:r2)'*ones(1,Nf)-patterns).^2.*(params(r2+[1:r2])'*ones(1,Nf))));
case 'Q2'
   %f = -sum(h.*(log(eps+(1./(1 + exp((params(1:r2) * patterns))))) + ...
   %   			 (1-h).*log(eps+1-(1./(1 + exp((params(1:r2)* patterns)))))));
   %f = -sum(h.*log(eps+exp(-0.5*sum((params(1:r2)'*ones(1,Nf)-patterns).^2.*(params(r2+[1:r2])'*ones(1,Nf))))) + ...
   %   (1-h).*log(eps+1-exp(-0.5*sum((params(1:r2)'*ones(1,Nf)-patterns).^2.*(params(r2+[1:r2])'*ones(1,Nf))))));
   f = LocBoostFunctions(params, 'gamma_kernel', patterns, [], [], params2);
   f = -sum(h.*log(eps+f) + (1-h).*log(eps+1-f));
case 'NewTestSet'
    %This section is used for labeling new data (especially of dimension > 2)
    %In this case, params is phi and params2 is theta
    phi                 = params;
    theta               = params2;
    [Dims, Nf]          = size(patterns);
    targets             = ones(1,Nf);
    patterns(Dims+1,:)  = ones(1,length(targets));
    
    Pdecision           = 0.5 * ones(1,Nf);
    %Pdecision           = LocBoostFunctions(theta(1,:), 'class_kernel', patterns, targets);
    
    for t = 2:size(params,1),
        Dgamma      = real(LocBoostFunctions(phi(t,:), 'gamma_kernel', patterns(1:Dims,:),[],[],Dims));  
        Dclass	    = LocBoostFunctions(theta(t,:), 'class_kernel', patterns, targets);
        Pdecision   = (1-Dgamma).*Pdecision + Dgamma.*Dclass;
    end
    
    f = Pdecision;
otherwise
   error ('Function type not recognized');
end


