%Code written by Harel Avissar and Danny Bickson
%Suplamentary material for the paper: 
%Distributed fault identification via non-parametric belief propagation.
%By D. Bickson, H. Avissar, D. Dolev, S. P. Boyd, A. Ihler and D. Baron.
%Submitted for publication. Aug 2009.
%Code available from http://www.cs.huji.ac.il/labs/danss/p2p/gabp/

%generates figure 5 in the above paper
function [res] = Generate_3D()
clear
avgT = 0;
totalR = 0;
LOFlag = 1;
K = 10;
n = 50;%15;
m = 30;%10;
Slevels = 2; %sparcity levels to check
Elevels = 1; %error levels to check
Sstart = 2;
Estart = 1;
algs = 8; %number of algorithems compare
% params are : 1- time      2- success rate         3- iterations to
% converge      4- number of times when found MAP solution better than real
% solution      5- ranking of solution          6- number of LO rounds
% 7- number of exact solution           8- number of times on top
INBP = 1; IIP = 2; ICoSaMP = 3; IGPSR = 4; IhIO = 5; ICS = 6; INON = 7; ILDLC = 8;%index for each alg
solutions = zeros(n, algs); %at each round
iterations = zeros(1, algs); %at each round
times = zeros(1, algs); %at each round
errors = zeros(1, algs); %for all rounds,will be shown every end of loop
%rand('state',76);%plot - 76
rand('state',1);
randn('state',1);

for er=Estart:Elevels
    
for sp=Sstart:Slevels
    
p = 0.05*er;%0.1
q = sp/10;%0.1;
skipped = 0;
test =1;
%tic;
for state=1:test
%load('H_n_121_d_3_1.mat','H');
    %%%H = rand(m,n);
    H = rand(m,n);
    %Hm = randn(m,n)*0.25+1;
    H(H<1-q) = 0;
    H(H>1-q+q/2) = -1;
    H(H>1-q) = 1;
    %H = H.*Hm;
    A = 2*H;%0.5*H;

    x = rand(n,1); % transmit a random vector of {-1,1}
    x(x>1-p) = 1;
    x(x<=1-p) = -1;
    [sortP, indP] = sort(x, 'descend');
    X = (x+1)/2;
    %x(x<0.5) = -1;
    %x(x>=0.5) = 1;
    sigma = 1;%0.4;
    Noise = randn(m,1)*sigma;

    y = H*x+Noise;
    Y = y+0.5*A*ones(n,1);
    
    prn = x;
    str = sprintf('Five biggest entries of x:[%d %d %d %d %d]', indP(1:5));
    disp(str);
    fprintf('Their values and corresponding convergence:\n');
    str = sprintf('[ %7.2f %7.2f %7.2f %7.2f %7.2f]', prn(indP(1:5)));
    disp(str);
    fprintf('NBP:\n');
    
    %% NBP
    max_rounds = 12;%12; % max number of iterations
    epsilon = 1e-20; % 
    boundx=n*(1.2*q+2*p);%n*q*2;%(1.2*q+2*p);
    model_order = 243;
    xx=1:model_order;  xx=xx-(model_order+1)/2;  xx=xx/max(xx);
    xx=xx*boundx; % values over which pdf is sampled
    pdf_prior= p*normpdf(xx,1,.1) + (1-p)*normpdf(xx,-1,.1);
    pdf_prior=pdf_prior/sum(pdf_prior);
    noise_prior=normpdf(xx,0,1); %noise prior
    noise_prior=noise_prior/sum(noise_prior);
    try
    tic;
    [xrecon, solutions(:,INBP), srecon, iterations(INBP),conNBP]=NBP(H,x',y,sigma, max_rounds, indP(1:5) ,epsilon,pdf_prior,noise_prior,xx,1,.1);
    times(INBP) = toc;
    iterations(INBP) = iterations(INBP)-1;
    catch ME
        skipped = skipped+1;
        errors(INBP) = errors(INBP)+1;
        %rethrow(ME);
        continue
    end
    fprintf('CS:\n');
    %% NON
    solutions(:,ICS) = solutions(:,INBP)*0;
    times(ICS) = 0;
    iterations(ICS) = 0;
    %% NBP_CS
    boundx=n*(1.2*q+2*p);
    model_order = 243;
    xx=1:model_order;  xx=xx-(model_order+1)/2;  xx=xx/max(xx);
    xx=xx*boundx; % values over which pdf is sampled
    pdf_prior= p*normpdf(xx,-1,3) + (1-p)*normpdf(xx,-1,.1);
    pdf_prior=pdf_prior/sum(pdf_prior);
    noise_prior=normpdf(xx,0,1); %noise prior
    noise_prior=noise_prior/sum(noise_prior);
    try
    tic;
    [xrecon, solutions(:,ICS), srecon, iterations(ICS),conCS]=NBP(H,x',y,sigma, max_rounds, indP(1:5) ,epsilon,pdf_prior,noise_prior,xx,1,.1);
    times(ICS) = toc;
    iterations(ICS) = iterations(ICS)-1;
    catch ME
        %skipped = skipped+1;
        solutions(:,ICS) = solutions(:,INBP)*0;
        times(ICS) = 0;
        iterations(ICS) = 0;
        errors(ICS) = errors(ICS)+1;
        %continue
    end
    fprintf('LDLC:\n');
    %% NBP_LDLC
    boundx=n*(1.2*q+2*p);
    model_order = 243;
    xx=1:model_order;  xx=xx-(model_order+1)/2;  xx=xx/max(xx);
    xx=xx*boundx; % values over which pdf is sampled
    pdf_prior= 0.5*normpdf(xx,1,.1) + 0.5*normpdf(xx,-1,.1);
    pdf_prior=pdf_prior/sum(pdf_prior);
    noise_prior=normpdf(xx,0,1); %noise prior
    noise_prior=noise_prior/sum(noise_prior);
    try
    tic;
    [xrecon, solutions(:,ILDLC), srecon, iterations(ILDLC),conLDLC]=NBP(H,x',y,sigma, max_rounds, indP(1:5) ,epsilon,pdf_prior,noise_prior,xx,1,.1);
    times(ILDLC) = toc;
    iterations(ILDLC) = iterations(ILDLC)-1;
    catch ME
        %skipped = skipped+1;
        solutions(:,ILDLC) = solutions(:,INBP)*0;
        times(ILDLC) = 0;
        iterations(ILDLC) = 0;
        errors(ILDLC) = errors(ILDLC)+1;
        %continue
    end
    fprintf('IP:\n');
    %% ----------- ARMAP solution only ------------------------------------
    try
    tic;
    [sol,iterations(IIP),conIP] = armap_gabp(A,X,Y,p*ones(n,1),sigma,1/(2*n),0,'verbose',0);
    solutions(:,IIP) = 2*(sol-0.5);
    times(IIP) = toc;
    catch ME
        skipped = skipped+1;
        errors(IIP) = errors(IIP)+1;
        rethrow(ME);
        %continue
    end
    for i=1:iterations(IIP)
        fprintf('[ %7.2f %7.2f %7.2f %7.2f %7.2f]\n', conIP(indP(1:5),i));
    end
    fprintf('\n');
    %% ----------- CoSaMP ------------------------------------
    try
    tic;
    [sol,iterations(ICoSaMP), conCoSaMP] = CoSaMP(Y,A,p*n,max_rounds,0);
    solutions(:,ICoSaMP) = 2*(sol-0.5);
    times(ICoSaMP) = toc;
    catch ME
        skipped = skipped+1;
        errors(ICoSaMP) = errors(ICoSaMP)+1;
        %rethrow(ME);
        continue
    end
    %% ----------- GPSR ------------------------------------
    try
    tic;
    [sol,iterations(IGPSR), conGPSR] = GPSR_BB(Y,A,0.3,0);
    solutions(:,IGPSR) = 2*(sol-0.5);
    times(IGPSR) = toc;
    catch ME
        skipped = skipped+1;
        errors(IGPSR) = errors(IGPSR)+1;
        %rethrow(ME);
        continue
    end
    %% ----------- hardIO ------------------------------------
    try
    [sol,iterations(IhIO), conhIO] = hard_l0_Mterm(Y,A,n,n*p,0);
    solutions(:,IhIO) = 2*(sol-0.5);
    times(IhIO) = toc;
    catch ME
        skipped = skipped+1;
        errors(IhIO) = errors(IhIO)+1;
        %rethrow(ME);
        continue
    end

end
end
end

plot3([1], [1], [1], 'O','MarkerSize', 8, 'MarkerEdgeColor', 'k','DisplayName','Solution','LineWidth',2);%, 'MarkerFaceColor', 'k');
hold on;
grid on;
plot3(conNBP(indP(1),:), conNBP(indP(2),:), conNBP(indP(3),:), '-xb','LineWidth',2,'DisplayName','NBP');
plot3(conIP(indP(1),:), conIP(indP(2),:), conIP(indP(3),:), '-xr','LineWidth',2,'DisplayName','IP');
plot3(conCS(indP(1),:), conCS(indP(2),:), conCS(indP(3),:), '-xy','LineWidth',2,'DisplayName','CS');
plot3(conLDLC(indP(1),:), conLDLC(indP(2),:), conLDLC(indP(3),:), '-xg','LineWidth',2,'DisplayName','LDLC');
legend('NBP','IP','CSBP','LDLC');

function [Xamb Lamb rounds]=Round_and_Local(s,K,H,p,n,y,sigma)
    %% Rounding
    [x_sort,ind_x] = sort(s,'descend');
    x_cand = []; l_cand = [];
    %x_sort
    for i = 1:min(max(K, 3*n*p),n)
        x_cur = -ones(n,1);
        x_cur(ind_x(1:i)) = 1;
        l_cur = L(H, p, n, x_cur, y, sigma);
        x_cand = [x_cand x_cur]; l_cand = [l_cand l_cur]; 
    end
    [l_sort,ind_l] = sort(l_cand,'ascend');
    X_amb = x_cand(:,ind_l(1:K)) ;%get ambiguity set
    L_amb = l_sort(1:K);
    
    iter = 0;
    if LOFlag
        %%LocalOpt
        x_cand = X_amb;
       l_cand = L_amb;
        EXIT_FLAG = 0;
        while(~EXIT_FLAG)
            x_cur = X_amb(:,1);x_best = x_cur;
            iter = iter+1;
            for j = 1:n
                x_cur(j) = -x_cur(j);
                if ismember(x_cur', x_cand', 'rows')
                    x_cur(j) = -x_cur(j);
                else
                    x_cand = [x_cand x_cur];
                    l_cand = [l_cand L(H, p, n, x_cur, y, sigma)];
                    x_cur(j) = -x_cur(j);
                end
            end
            [l_sort,ind_l] = sort(l_cand,'ascend');
            X_amb = x_cand(:,ind_l(1:K)); %get ambiguity set
            L_amb = l_sort(1:K);
            if all(x_best==X_amb(:,1))
                EXIT_FLAG = 1;
            end
        end
    end
    Xamb = X_amb;
    Lamb = L_amb;
    rounds = iter;    
end

function Likelihood=L(H, p, n, x, y, sigma)
    Likelihood = (1/(2*sigma^2))*trace(H'*H*x*x')+(log((1-p)/p)*ones(n,1)-(1/sigma^2)*H'*y)'*x;
end

%fprintf('GaBP success=%6.24  \n',sum(x == sign(inv(H)*y))/121);
end