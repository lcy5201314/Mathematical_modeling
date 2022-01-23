

% script for creating figure 2.
% Algorithm is described in the paper:
% "Fixing the convergence of the GaBP algorithm" by
% J. K. Johnson, D. Bickson and D. Dolev
% In ISIT 2009
% http://arxiv.org/abs/0901.4192
% Code written by Danny Bickson
% December 2008. Updated May 2009.
% 
% input: gamma - ratio of diagonal increase
function [] = multiuser_detection_fix2(gamma)
rand('state',1);
randn('state',1);

addpath('..');

if (nargin < 1)
    gamma = 1;
end

N=256; %number of chips per symbol
K=64; %number of active users.
A=eye(K); %assuming equal power for now.

sigma2=1; % from dB=10*log10(1/sigma2)
b= rand(K,1); %transmitted bits
b(b < 0.5) = -1;
b(b >= 0.5) = 1;

S = rand(N,K); %random spreading CDMA 
S(S < 0.5) = -1;
S(S >= 0.5) = 1;
S = S/sqrt(N);

a11 = eig(eye(K) - abs(S'*S));
disp(['spectral radius is: ' , num2str(max(abs(a11)))]);

n=sqrt(sigma2)*randn(N,1); %AWGN
r=S*A*b+n; %received signal
y=S'*r; %matched filter output
        
M=S'*S + sigma2*eye(K);% correlation matrix
b_est2 = inv(M)*y;
inv_err=sum(b~=(b_est2 > 0))/K;
disp(['err using direct inversion is: ', num2str(inv_err)]);

epsilon_inner_loop = 1e-6;
epsilon_outer_loop = 1e-3;

MAX_ROUNDS = 1000;
xj = zeros(K,MAX_ROUNDS);
xj(:,1) = y;
xj2 = zeros(K,MAX_ROUNDS);
xj2(:,1) = y;

diagonal_loading = max(sum(abs(M)) - max(diag(M))) * gamma;
Minc = M + eye(length(M)) * diagonal_loading;
Minc2 = M + eye(length(M)) * diagonal_loading/5;


a11 = eig(eye(length(M)) - Minc*inv(diag(diag(Minc))));
a22 = eig(eye(length(M)) - Minc2*inv(diag(diag(Minc2))));

norms = ones(1,MAX_ROUNDS);
norms2 = ones(1,MAX_ROUNDS);
for l=2:MAX_ROUNDS
    % This is a single Newton step
    % Algorithm is described in
    % Linear Detection via Belief Propagation. Danny Bickson, Danny Dolev, Ori
    % Shental, Paul H. Siegel and Jack K. Wolf. In the 45th Annual Allerton Conference on Communication, Control, 
    % and Computing, Allerton House, Illinois, Sept. 07.

    [b_est,J,r1] = asynch_GBP(Minc, y - M*xj(:,l-1), MAX_ROUNDS, epsilon_inner_loop);
    xj(:,l) = xj(:,l-1) + b_est';
    e = norm(y - M*xj(:,l));
    disp(['error norm 1 for round ', num2str(l), ' is ', num2str(e)]);

    [b_est2,J2,r12] = asynch_GBP(Minc2, y - M*xj2(:,l-1), MAX_ROUNDS, epsilon_inner_loop);
    xj2(:,l) = xj2(:,l-1) + b_est2';
    e2 = norm(y - M*xj2(:,l));
    disp(['error norm 2 for round ', num2str(l), ' is ', num2str(e2)]);

    if (e < epsilon_outer_loop)
        break;
    end
    
    norms(l) = e;
    norms2(l) = e2;
end

       
figure;
hold on;
title('Convergence of fixed GaBP iteration with n=256,k=64','FontSize',16);
xlabel('Newton step','FontSize',14);
ylabel('Error norm','FontSize',14);
legend('\rho = ', '\rho = ');
semilogy(3:l-1, norms(3:l-1),'b*');
semilogy(3:l-1, norms2(3:l-1),'g+');

disp(['diagonal loading is ' num2str( diagonal_loading) ' gamma = ' num2str(gamma)]);
disp(['spectral radius is: ' , num2str(max(abs(a11)))]);
disp(['spectral radius 2 is: ' , num2str(max(abs(a22)))]);
end
