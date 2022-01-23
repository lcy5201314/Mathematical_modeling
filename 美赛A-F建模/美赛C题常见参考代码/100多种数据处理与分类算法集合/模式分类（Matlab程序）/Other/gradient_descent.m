function Min = gradient_descent(a, theta, eta, fun)

%	Minimize a function using the basic gradient descent algorithm
%
%	Inputs:
%		a		- Initial search point
%		theta   - Convergence criterion
%		eta	    - Descent rate
%		fun	    - The function to minimize
%
%	Outputs:
%		Min	    - The minimum point

gradJ			= 1e3;
step_size	= 1e-2;
D				= length(a);

while (gradJ > theta),
   %Compute the gradient of J in each direction
   gradJi = zeros(1,D);
   for i = 1:D,
      delta         = zeros(size(a));
      delta(i)      = step_size;
      y_minus		= feval(fun, a - delta);
      y_plus		= feval(fun, a + delta);
      gradJi(i)     = (y_plus - y_minus)/(2*step_size);
   end
   
   %a <- a - eta*gradJ
   a		= a - eta*gradJi;
   
   gradJ = sum(abs(gradJi))*eta;
end

Min = a;