%Code written by Danny Bickson
%Code is found on: http://www.cs.huji.ac.il/labs/danss/p2p/gabp/index.html
%Example of solving boolean least squares via NBP
%The problem is argmin_x ( norm (Hx - y))  s.t. x \in {-1,1}
%Ulike ldlc G is the sparse ENCODING matrix (and not the decoding)

%function [sigmas,success] = loop_LDLC_example()
clear

%load('H_N_121_d_3_1.mat','H');
%load('H_N_961_d_7.mat','H');
%sum(H)
%H(:,1)
%H(H>0) = 1;
%H(H<0) = -1;
% for i=1:1000
%     rand('state',i);
%     R = rand(961,961);
%     R(R>0.4) = 1;
%     R(R<0.6) = 0;
%     H = H.*R;
%     %sum(H)
%     G = H;
%     n=size(H,1);
%     detG=det(G);
%     if (degG~=0)
%         i
%     end
% end

n = 100; % size of the matrix
d = 3; % number of non zeros

rand('state',361);
H = zeros(n,n);
for j=1:n
    R = rand(2,d);
    for k=1:d
        while H(j,ceil(R(2,k)*n))~=0
            R = rand(2,d);
        end
        if R(1,k)<0.5
            H(j,ceil(R(2,k)*n))=1;
        else
            H(j,ceil(R(2,k)*n))=-1;
        end
    end
end
G = H;
n=size(H,1);
detG=det(G)
max_sigma_squared = nthroot(detG^2, n)/(2*pi*exp(1));
disp(['max sigma^2 for channel capacity is ', num2str(max_sigma_squared)]);

retry=10;
max_instance = 250;
max_rounds = 10; % max number of iterations

sigmas=linspace(max_sigma_squared/5, max_sigma_squared, retry);
success = zeros(1, retry);

epsilon = 1e-30; % 
boundx=10;

%---------------
% compute signal node prior
%---------------
model_order = 729;%1543;
xx=1:model_order;  xx=xx-(model_order+1)/2;  xx=xx/max(xx);
xx=xx*boundx; % values over which pdf is sampled

pdf_prior= normpdf(xx,-1,.2) + normpdf(xx,1,.2);
%pdf_prior= normpdf(xx,0,1);
pdf_prior=pdf_prior/sum(pdf_prior);

noise_prior=normpdf(xx,0,1); %noise prior
noise_prior=noise_prior/sum(noise_prior);

for loop_count=1:retry
    rand('state',loop_count);
    randn('state',loop_count);
    total = 0;
    for instance=1:max_instance
       
       x = rand(n,1); % transmit a random vector of {-1,1}
        x(x<0.5) = -1;
        x(x>=0.5) = 1;
        sigma = sigmas(loop_count);

        y = G*x+randn(n,1)*sqrt(sigma);

        disp (['STARTING THE DECODER... with sigma ', num2str(sigma), ' distance from capcity ' , num2str(-10*log10(sigma/max_sigma_squared)), ' dB']);
         [xrecon, mrecon, srecon]=NBP(G,x',y,...
             sigma, max_rounds, [ 1 10 21 24] ,...
             epsilon,pdf_prior,noise_prior,xx);


        fprintf('[   %d   %1d   %d   %d]\n', x([1 10 21 24]));
        sucrate = sum((((mrecon>0)*2)-1)==x')/n;
        fprintf('success=%6.4f \n',sucrate);
        total = total + sucrate;
    end
    total = total/max_instance;
    disp(['Avg performance with ', num2str(sigma),' is ', num2str(total)]);
    success(loop_count) = total;
end
save ret.mat sigmas success;
%end