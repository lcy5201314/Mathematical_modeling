

% Function that demonstrates how to fix the convergence of the GaBP
% algorithm, in the multiuser detection secnario

% Algorithm is described in the paper:
% "Fixing the convergence of the GaBP algorithm" by
% J. K. Johnson, M. Chetrkov, D. Bickson and D. Dolev
% Submitted to ISIT 2009
% Code written by Danny Bickson
% December 2008
function [x] = fix_GaBP(A,b,max_iter,epsilon)

dloading = max(sum(abs(A)) - max(diag(A)));
Minc = A + eye(length(A)) *  dloading;
disp(['diagonal loading is ' num2str(dloading)]);
if (dloading == 0)
   warning('No need to use double loop construction, since GaBP converges anyway. aborting');
   return;
end

xj = b;

for l=2:max_iter
    % This is a single Newton step
    % Algorithm is described in
    % Linear Detection via Belief Propagation. Danny Bickson, Danny Dolev, Ori
    % Shental, Paul H. Siegel and Jack K. Wolf. In the 45th Annual Allerton Conference on Communication, Control, 
    % and Computing, Allerton House, Illinois, Sept. 07.

   [xj,J,r1] = asynch_GBP(Minc, b + dloading*xj, max_iter, epsilon);
    xj = xj';
    disp(['error norm for round ', num2str(l), ' is ', num2str(norm(b - A*xj))]);
end

       
x = xj;


end
