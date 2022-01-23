function [theta, theta_debias, objective, times, debias_start, mses, num_nz_theta, Num_iteration, psnr_val] = My_GPSR_BB(y, A, tau, varargin)
% This function provides a solver the convex problem 
% arg min_theta = 0.5*|| y - A theta ||_2^2 + tau ||theta||_1
% %  ===== Required inputs =============
%  y: 1D vector or 2D array (image) of observations
%    
%  A: if y and theta are both 1D vectors, A can be a 
%     k*n (where k is the size of y and n the size of theta)
%     matrix or a handle to a function that computes
%     products of the form A*v, for some vector v.
%     In any other case (if y and/or theta are 2D arrays), 
%     A has to be passed as a handle to a function which computes 
%     products of the form A*x; another handle to a function 
%     AT which computes products of the form A'*x is also required 
%     in this case. The size of theta is determined as the size
%     of the result of applying AT.
%
%  tau: usually, a non-negative real parameter of the objective 
%       function (see above). It can also be an array, the same 
%       size as theta, with non-negative entries; in this case,
%       the objective function weights differently each element 
%       of theta, that is, it becomes
%       0.5*|| y - A theta ||_2^2 + tau^T * abs(theta)
%
%  ===== Optional inputs =============
%  'AT'    = function handle for the function that implements
%            the multiplication by the conjugate of A, when A
%            is a function handle. If A is an array, AT is ignored.
%
%  'StopCriterion' = type of stopping criterion to use
%                    0 = algorithm stops when the relative 
%                        change in the number of non-zero 
%                        components of the estimate falls 
%                        below 'ToleranceA'
%        any non-zero = algorithm stops when the relative 
%                       change in the objective function 
%                       falls below 'ToleranceA', for the
%                       first phase and 'ToleranceD' for 
%                       debiasing phase.
%                       Default = 0
%
%  'ToleranceA' = stopping threshold; Default = 0.01
% 
%  'Debias'     = debiasing option: 1 = yes, 0 = no.
%                 Default = 0.
%
%  'ToleranceD' = stopping threshold for the debiasing phase:
%                 Default = 0.0001.
%                 If no debiasing takes place, this parameter,
%                 if present, is ignored.
%
%  'MaxiterA' = maximum number of iterations allowed in the
%               main phase of the algorithm.
%               Default = 200
%
%  'MiniterA' = minimum number of iterations performed in the
%               main phase of the algorithm.
%               Default = 5
%
%  'MaxiterD' = maximum number of iterations allowed in the
%               debising phase of the algorithm.
%               Default = 200
%
%  'MiniterD' = minimum number of iterations to perform in the
%               debiasing phase of the algorithm.
%               Default = 5
%
%  'Initialization' must be one of {0,1,2,array}
%               0 -> Initialization at zero. 
%               1 -> Random initialization.
%               2 -> initialization with A'*y.
%           array -> initialization provided by the user.
%               Default = 0;
%             
%  'TrueTheta' = if the true underlying theta is passed in 
%                this argument, MSE plots are generated.
%
%  'AlphaMin' = the alphamin parameter of the BB method.
%               (See the paper for details)
%               Default = 1e-30;
%
%  'AlphaMax' = the alphamax parameter of the BB method.
%               (See the paper for details)
%               Default = 1e30;
%
% ============ Outputs ==============================
%   theta = solution of the main algorithm
%
%   theta_debias = solution after the debiasing phase;
%                  if no debiasing phase took place, this
%                  variable is empty, theta_debias = [].
%
%   objective = sequence of values of the objective function
%
%   times = CPU time after each iteration
%
%   debias_start = iteration number at which the debiasing 
%                  phase started. If no debiasing took place,
%                  this variable is returned as zero.
%
%   mses = sequence of MSE values, with respect to TrueTheta,
%          if it was given; if it was not given, mses is empty,
%          mses = [].
% ========================================================
% test for number of required parametres
if (nargin-length(varargin)) ~= 3
     error('Wrong number of required parameters');
end

% Set the defaults for the optional parameters
stop_objective = 0;
tolA = 0.01;
tolD = 0.0001;
debias = 0;
maxiter = 200;
maxiter_debias = 200;
miniter = 0;
miniter_debias = 5;
init = 0;
alphamin = 1e-30;
alphamax = 1e30;
compute_mse = 0;
AT = 0;

w1 = 0.9;
w2 = 0.1;

% Set the defaults for outputs that may not be computed
debias_start = 0;
theta_debias = [];
mses = [];

% Read the optional parameters
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    switch varargin{i}
     case 'StopCriterion'
       stop_objective = varargin{i+1};
     case 'ToleranceA'       
       tolA = varargin{i+1};
     case 'ToleranceD'
       tolD = varargin{i+1};
     case 'Debias'
       debias = varargin{i+1};
     case 'MaxiterA'
       maxiter = varargin{i+1};
     case 'MaxiterD'
       maxiter_debias = varargin{i+1};
     case 'MiniterA'
       miniter = varargin{i+1};
     case 'MiniterD'
       miniter_debias = varargin{i+1};
     case 'Initialization'
       if prod(size(varargin{i+1})) > 1   % we have an initial theta
          init = 33333;    % some flag to be used below
          theta = varargin{i+1};
       else 
          init = varargin{i+1};
       end
     case 'TrueTheta'
       compute_mse = 1;
       true = varargin{i+1};
     case 'AlphaMin'
       alphamin = varargin{i+1};
     case 'AlphaMax'
       alphamax = varargin{i+1};
     case 'AT'
       AT = varargin{i+1};
     case 'SI'
       SI = varargin{i+1};
     case 'OriX'
       x0 = varargin{i+1};
       [m n] = size(x0);
       [h0, h1, l0, l1] = filter9_7();
       L = floor(log2(m)) - 3;    
       x = x0(:);
       xmean = mean(x);
     case 'OriSI'
       s0 = varargin{i+1};
     case 'Num_of_Measurements'
       K = varargin{i+1};
       CS_Rate = K / (m*n);
       %fprintf(1, 'CS rate = %.2f\n', CS_Rate);
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized option: ''' varargin{i} '''']);
    end;
  end;
end
%%%%%%%%%%%%%%

% if A is a function handle, we have to check presence of AT,
if isa(A, 'function_handle') & ~isa(AT,'function_handle')
   error(['The function handle for transpose of A is missing']);
end 

% if A is a matrix, we find out dimensions of y and theta,
% and create function handles for multiplication by A and A',
% so that the code below doesn't have to distinguish between
% the handle/not-handle cases
if ~isa(A, 'function_handle')
   AT = @(x) A'*x;
   A = @(x) A*x;
end
% from this point down, A and AT are always function handles.

% start the clock
t0 = cputime;

% Precompute A'*y since it'll be used a lot
Aty = AT(y);

% Initialization
switch init
    case 0   % initialize at zero, using AT to find the size of theta
       theta = AT(zeros(size(y)));
    case 1   % initialize randomly, using AT to find the size of theta
       theta = randn(size(AT(zeros(size(y)))));
    case 2   % initialize theta0 = A'*y
       theta = Aty; 
    case 33333
       % initial theta was given as a function argument; just check size
       if size(A(theta)) ~= size(y)
          error(['Size of initial theta is not compatible with A']); 
       end
    otherwise
       error(['Unknown ''Initialization'' option']);
end

% now check if tau is an array; if it is, it has to 
% have the same size as theta
if prod(size(tau)) > 1
   try,
      dummy = theta.*tau;
   catch,
      error(['Parameter tau has wrong dimensions; it should be scalar or size(theta)']),
   end
end

% if the true theta was given, check its size
if compute_mse & (size(true) ~= size(theta))  
   error(['Initial theta has incompatible size']); 
end

% initialize u, and v
u =  theta.*(theta >= 0);
v = -theta.*(theta <  0);

% define the indicator vector or matrix of nonzeros in theta
nz_theta = (theta ~= 0.0);

% Compute and store initial value of the objective function
resid =  y - A(theta);
f11 = 0.5*(resid(:)'*resid(:));
f12 = sum(tau(:).*u(:)) + sum(tau(:).*v(:));
prev_f =  f11 + f12;         

%%%%%
    f1 = prev_f;
    prev_f1 = f1;
    f2 = (sum(sum((theta - SI).^2)));
    prev_f = w1*f1 + w2*f2;

   s_theta = reshape(theta, m, n);
   xr = idwt2d(s_theta, l0, l1, L);
   xr = xr + xmean;
   psnr_val  = psnr(x0, xr);
   %Lap_x0_xr = sqrt(2) / std(x0(:) - xr(:));
   Lap_s0_xr = sqrt(2) / std(s0(:) - xr(:));
   %Lap_x0_s0 = sqrt(2) / std(x0(:) - s0(:));
   %fprintf(1,'Iter=  -1, obj = %.4f, (f11, f12, f2) = (%.2f, %.2f, %.2f), PSNR = %.2f\n', prev_f, f11, f12, f2, psnr_val);
   %fprintf(1, '(Lap_x0_s0, Lap_s0_xr, Lap_x0_xr) = (%.2f, %.2f, %.2f)\n', Lap_x0_s0, Lap_s0_xr, Lap_x0_xr);
%%%%%

times(1) = cputime - t0;
objective(1) = prev_f;

if compute_mse
   mses(1) = sum(sum((theta-true).^2));
end

% Compute the initial gradient and the useful 
% quantity resid_base
resid_base = y - resid;
temp = AT(resid_base);
term = temp - Aty;
gradu =   term + tau;
gradv =  -term + tau;

% projection and computation of search direction vector
du = max(u - gradu, 0.0) - u;
dv = max(v - gradv, 0.0) - v;

% initialization of alpha
alpha = 1/max(max(abs(du(:))),max(abs(dv(:))));
alphas(1) = alpha;

% control variable for the outer loop and iteration counter
% cont_outer = (norm(projected_gradient) > 1.e-5) ;

cont_outer = 1;
iter = 0;

%%%%%
if prev_f == 0
    f = prev_f;
    theta = SI;
    cont_outer = 0;
end
%%%%%

while cont_outer

  theta_previous = theta;
  psnr_previous = psnr_val;
  previous_Lap_s0_xr = Lap_s0_xr;
  
  % compute gradient
  temp = AT(resid_base);
  
  term  =  temp - Aty;
  gradu =  term + tau;
  gradv = -term + tau;
  
  % projection and computation of search direction vector
  du = max(u - alpha*gradu, 0.0) - u;
  dv = max(v - alpha*gradv, 0.0) - v;
  dtheta = du-dv;
  old_u = u;
  old_v = v;
  
  % calculate useful matrix-vector product involving dtheta
  auv = A(dtheta);
  dGd = auv(:)'*auv(:);
  if (dGd > 0)
    lambda0 = - (gradu(:)'*du(:) + gradv(:)'*dv(:))/dGd;
    if (lambda0 <= 0)
      % if we get to here, there's something wrong
      alpha = alpha / 2; 
      fprintf(1,'alpha = %g\n',alpha);
      lambda = 0;
    else
     lambda = min(lambda0,1);
     u = old_u + lambda * du;
     v = old_v + lambda * dv;
     uvmin = min(u,v);
     u = u - uvmin; 
     v = v - uvmin; 
     theta = u - v;
    
    % update the nonzero indicator vector
     nz_theta_prev = nz_theta;
     nz_theta = (theta~=0.0);
     num_nz_theta = sum(nz_theta(:));
     num_changes_active = sum(nz_theta(:)~=nz_theta_prev(:));

     resid = y - resid_base - lambda*auv;
     f11 = 0.5*(resid(:)'*resid(:));
     f12 = sum(tau(:).*u(:)) + sum(tau(:).*v(:));
     f = f11 + f12;
        
  %%%%%
      f1 = f;
      f2 = (sum(sum((theta- SI).^2)));
      f = w1*f1 + w2*f2;
  %%%%%   
    
     dd  = du(:)'*du(:) + dv(:)'*dv(:); 
     alpha = min(alphamax,max(alphamin,dd/dGd));
     resid_base = resid_base + lambda*auv;
    end  
    
  else %  (dGd < 0)
      
      % something wrong if we get to here
      fprintf(1,' dGd=%12.4e, something is wrong\n', dGd);
      alpha = alphamax;     
%       resid_base = resid_base + auv;
  end

   s_theta = reshape(theta, m, n);
   xr = idwt2d(s_theta, l0, l1, L);
   xr = xr + xmean;
   psnr_val  = psnr(x0, xr); 
   %Lap_x0_xr = sqrt(2) / std(x0(:) - xr(:));
   Lap_s0_xr = sqrt(2) / std(s0(:) - xr(:));
   %Lap_x0_s0 = sqrt(2) / std(x0(:) - s0(:));
   %fprintf(1,'Iter=%4d, obj = %.4f, (f11, f12, f2) = (%.2f, %.2f, %.2f), PSNR = %.2f\n', iter, f, f11, f12, f2, psnr_val);
   %fprintf(1, '(Lap_x0_s0, Lap_s0_xr, Lap_x0_xr) = (%.2f, %.2f, %.2f)\n', Lap_x0_s0, Lap_s0_xr, Lap_x0_xr);

% evaluate projected gradient quantity of Dai-Fletcher
%   projected_gradient = [min(gradu,0.0).*(u == 0.0) + gradu.*(u ~= 0.0);
%     min(gradv,0.0).*(v==0.0) + gradv.*(v~=0.0)];
  
  if ~stop_objective
     diff_f = f - prev_f; 
     ObjChangeRate = abs(diff_f) / prev_f;
     LapChangeRate = abs(Lap_s0_xr - previous_Lap_s0_xr) / previous_Lap_s0_xr;
     %fprintf(1, '(ObjChangeRate, LapChangeRat) = (%.2f, %.2f)\n', ObjChangeRate, LapChangeRate);
     
     if CS_Rate <= 0.2
        if (LapChangeRate >= 0.9) & (iter > miniter)
            cont_outer = 0;
        end
     elseif (CS_Rate > 0.2) & (CS_Rate <= 0.7)
        if (LapChangeRate <= 0.05) & (iter > miniter)
            cont_outer = 0;
            if ObjChangeRate >= 0.01
               theta = theta_previous;
               psnr_val = psnr_previous;                  
            end
        elseif (diff_f >= 0) & (iter > miniter)
            cont_outer = 0;
            if ObjChangeRate >= 0.01
               theta = theta_previous;
               psnr_val = psnr_previous;                  
            end
        end
     elseif (CS_Rate > 0.7)
         if ObjChangeRate <= 0.001
             cont_outer = 0;
             theta = theta_previous;
             psnr_val = psnr_previous;   
         end
     end
     if iter > maxiter
        cont_outer = 0;
     end
  else
     crit = abs(f-prev_f)/(prev_f*tolA);
     %fprintf(1,'Iter=%4d, obj = %.4f, (f1, f2) = (%.2f, %.2f), PSNR = %.2f\n', iter, f, f1, f2, psnr_val);
  end

  iter = iter + 1;
  prev_f = f;
  prev_f1 = f1;
  objective(iter) = f;
  
  times(iter) = cputime - t0;
  alphas(iter) = alpha;
  
  if compute_mse
     err = true - theta;
     mses(iter) = (err(:)'*err(:));
  end
  
   if ~stop_objective
       %cont_outer = ((iter <= maxiter) & ...
       %             ((num_changes_active > tolA*num_nz_theta)|...
       %             (lambda<1.0)));
   else 
       cont_outer = ((iter <= maxiter) & ...
                    ((crit > 1)|(lambda<1.0)));
   end
  
  %if psnr_val - psnr_previous < 0
     %theta = theta_previous;
     %psnr_val = psnr_previous;
     %cont_outer = 0;
  %end
end % end of the main loop of the GPBB algorithm
Num_iteration = iter + 2;
% Printout results
%fprintf(1,'\nFinished the main algorithm!\nResults:\n');
%fprintf(1,'||A theta - y ||_2 = %.2f\n',resid(:)'*resid(:));
%fprintf(1,'||theta||_1 = %.2f\n',sum(abs(theta(:))));
%fprintf(1,'Objective function = %.2f\n',f);
%fprintf(1,'Number of non-zero components = %d\n\n', num_nz_theta);
%fprintf(1, 'Final psnr_val = %.2f\n\n', psnr_val);
%pause;