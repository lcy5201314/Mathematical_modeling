

% Function that demonstrates how to fix the convergence of the GaBP
% algorithm, in the multiuser detection secnario

% Algorithm is described in the paper:
% "Fixing the convergence of the GaBP algorithm" by
% J. K. Johnson, D. Bickson and D. Dolev
% In ISIT 2009
% http://arxiv.org/abs/0901.4192
%
% Code written by Danny Bickson
% December 2008
% 
% input: gamma - ratio of diagonal increase
function [inner,outer,rho] = multiuser_detection_fix3(gamma)
rand('state',1);
randn('state',1);

addpath('..');

if (nargin < 1)
    gamma = 1;
end

N=256; %number of chips per symbol
K=128; %number of active users.
 %assuming equal power for now.

sigma2=1; % from dB=10*log10(1/sigma2)
b= rand(K,1); %transmitted bits
b(b < 0.5) = -1;
b(b >= 0.5) = 1;

S = rand(N,K); %random spreading CDMA 
S(S < 0.5) = -1;
S(S >= 0.5) = 1;
S = S/sqrt(N);

A = eye(K) - (sigma2*eye(K) + S'*S);
A = A^-1/2 * A * A^-1/2;
a11 = eig(abs(A));
rho = max(abs(a11));
disp(['spectral radius is: ' , num2str(rho)]);

n=sqrt(sigma2)*randn(N,1); %AWGN
r=S*b+n; %received signal
y=S'*r; %matched filter output
        
M=S'*S;% correlation matrix
b_est2 = inv(M)*y;
inv_err=sum(b~=(b_est2 > 0))/K;
disp(['err using direct inversion is: ', num2str(inv_err)]);

epsilon_inner_loop = 1e-6;
epsilon_outer_loop = 1e-3;

MAX_ROUNDS = 1000;
xj = zeros(K,MAX_ROUNDS);
xj(:,1) = y;

diagonal_loading = max(sum(abs(M)) - max(diag(M))) * gamma;
Minc = M + eye(K) * diagonal_loading;
Minc2 = eye(K) - M;
Minc2 = Minc2^-1/2*Minc2*Minc2^-1/2;
a11 = eig(abs(Minc2));
rho = max(abs(a11));
disp(['spectral radius of Minc is: ' , num2str(rho)]);

inner = 0;
norms = ones(1,MAX_ROUNDS);

for l=2:MAX_ROUNDS
    % This is a single Newton step
    % Algorithm is described in
    % Linear Detection via Belief Propagation. Danny Bickson, Danny Dolev, Ori
    % Shental, Paul H. Siegel and Jack K. Wolf. In the 45th Annual Allerton Conference on Communication, Control, 
    % and Computing, Allerton House, Illinois, Sept. 07.

    [b_est,J,temp] = asynch_GBP(Minc, y - M*xj(:,l-1), MAX_ROUNDS, epsilon_inner_loop);
    inner = inner + temp;
    xj(:,l) = xj(:,l-1) + b_est';
    
    e = norm(y - M*xj(:,l));
    disp(['error norm for round ', num2str(l), ' is ', num2str(e)]);
 
    norms(l) = e;
    %e = norm(y - M*xj(:,l));
    %disp(['error norm for round ', num2str(l), ' is ', num2str(e)]);

    if (e < epsilon_outer_loop)
		outer = l;
        inner = inner / l;
        break;
        
    end
    
    
   
end

       
% figure;
% hold on;
% semilogy(3:l-1, norms(3:l-1));
% title('Convergence of fixed GaBP iteration with n=256,k=64');
% xlabel('Newton step');
% ylabel('Error norm');
%legend('Without damping', 'With damping, gamma=0.5');

disp(['diagonal loading is ' num2str( diagonal_loading) ' gamma = ' num2str(gamma)]);
disp(['spectral radius is: ' , num2str(max(abs(a11)))]);
end
