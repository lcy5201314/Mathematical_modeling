function Min = Newton_descent(a, theta, fun)

%	Minimize a function using the Newton descent algorithm
%
%	Inputs:
%		a		- Initial search point
%		theta   - Convergence criterion
%		fun	    - The function to minimize
%
%	Outputs:
%		Min	    - The minimum point

gradJ			= 1e3;
step_size	= 1e-2;
D				= length(a);

while (gradJ > theta),
   %Compute the gradient of J in each direction, as well as the Hessian matrix
   gradJi = zeros(1,D);
   H		 = zeros(D,D);
   
   for i = 1:D,
      delta         = zeros(size(a));
      delta(i)      = step_size;
      y_minus		= feval(fun, a - delta);
      y_plus		= feval(fun, a + delta);
      gradJi(i)	    = (y_plus - y_minus)/(2*step_size);
   end
   
   for i = 1:D,
      for j = i:D,
          deltai    = zeros(size(a));
          deltai(i) = step_size;
          deltaj    = zeros(size(a));
          deltaj(j) = step_size;
          H(i,j)    = (feval(fun, a+deltai+deltaj) - feval(fun, a+deltai-deltaj) - ...
                       feval(fun, a-deltai+deltaj) + feval(fun, a-deltai-deltaj))/(4*step_size^2);
          H(i,j)    = H(j,i);
      end
   end
   
   update = inv(H)*gradJi';
   
   if isnan(sum(update)),
       break
   else
       %a <- a - inv(H)*gradJ
       a		= a - (inv(H)*gradJi')';

       gradJ = sum(abs(update));
   end
   
   disp(num2str(a))
end

Min = a;