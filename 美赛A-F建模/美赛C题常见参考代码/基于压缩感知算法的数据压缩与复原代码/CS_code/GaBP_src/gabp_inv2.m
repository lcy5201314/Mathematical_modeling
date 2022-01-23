%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.

%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.

%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implementation of the Gaussian BP algorithm, as given in: 
% Linear Detection via Belief Propagation
% By Danny Bickson, Danny Dolev, Ori Shental, Paul H. Siegel and Jack K. Wolf.
% In the 45th Annual Allerton Conference on Communication, Control and Computing, Allerton House, Illinois, 
% Sept. 07'
%
%
% Written by Danny Bickson.
% updated: 24-Dec-2008
%
% This code computes inv(A), given A.

% input: A - square matrix nxn
% max_iter - maximal number of iterations
% epsilon - convergence threshold
% output: Ainv - the matrix inv(A)
%
% note: this matlab code is NOT optimized, and serves as a demonstration
function [Ainv] = gabp_inv2(A, max_iter_in, max_iter_out, epsilon, quiet)

n = length(A);

for k=1:n
    b = zeros(n,1);
    b(k) = 1;
    x = fix_GaBP(A,b,max_iter_in,max_iter_out,epsilon,1,quiet);
    Ainv(:,k) = x;
end

end
