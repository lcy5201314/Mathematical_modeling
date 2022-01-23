%Code written by Harel Avissar and Danny Bickson
%Suplamentary material for the paper: 
%Distributed fault identification via non-parametric belief propagation.
%By D. Bickson, H. Avissar, D. Dolev, S. P. Boyd, A. Ihler and D. Baron.
%Submitted for publication. Aug 2009.
%Code available from http://www.cs.huji.ac.il/labs/danss/p2p/gabp/

% code for creating figure 3

function [res] = Generate_SuccRate()
clear
avgT = 0;
totalR = 0;
LOFlag = 0;
K = 10;
n = 100;%15;
m = 50;%10;
Slevels = 1; %TODO %sparsity levels to check
Elevels = 10; %error levels to check
Sstart = 1;
Estart = 1;
algs = 8; %number of algorithems compare
params = 8; %number of parameters to monitor on each alg
% params are : 1- time      2- success rate         3- iterations to
% converge      4- number of times when found MAP solution better than real
% solution      5- ranking of solution          6- number of LO rounds
% 7- number of exact solution           8- number of times on top
INBP = 1; IIP = 2; ICoSaMP = 3; IGPSR = 4; IhIO = 5; ICS = 6; INON = 7; ILDLC = 8;%index for each alg
names = ['NBP   ';'IP    ';'CoSaMP';'GPSR  ';'hardIO';'CS    ';'NON   ';'LDLC  '];
colors = ['b', 'r', 'c', 'm', 'k', 'g', 'y', 'w'];
results = zeros(Elevels,Slevels,algs,params); %of comparisons
solutions = zeros(n, algs); %at each round
best = zeros(1, algs); %at each round
iterations = zeros(1, algs); %at each round
times = zeros(1, algs); %at each round
errors = zeros(1, algs); %for all rounds,will be shown every end of loop
%rand('state',76);%plot - 76
rand('state',1111);
randn('state',1111);

for er=Estart:Elevels
    
for sp=Sstart:Slevels
    
p = 0.05*er;%0.1
q = sp/10;%0.1;
skipped = 0;
test =500; %%TODO
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
    X = (x+1)/2;
    %x(x<0.5) = -1;
    %x(x>=0.5) = 1;
    sigma = 1;%0.4;
    Noise = randn(m,1)*sigma;

    y = H*x+Noise;
    Y = y+0.5*A*ones(n,1);
    
    prn = x;
    str = sprintf('[ %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f]', prn([ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]));
    %%%disp(str);
    %%%fprintf('\n');
    
    
    %% NBP
    max_rounds = 8;%12; % max number of iterations
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
    [xrecon, solutions(:,INBP), srecon, iterations(INBP),conNBP]=NBP(H,x',y,sigma, max_rounds, [] ,epsilon,pdf_prior,noise_prior,xx,0,1);
    times(INBP) = toc;
    iterations(INBP) = iterations(INBP)-1;
    catch ME
        skipped = skipped+1;
        errors(INBP) = errors(INBP)+1;
        %rethrow(ME);
        continue
    end
    %% NON
    solutions(:,INON) = solutions(:,INBP)*0;
    times(INON) = 0;
    iterations(INON) = 0;
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
    [xrecon, solutions(:,ICS), srecon, iterations(ICS),conCS]=NBP(H,x',y,sigma, max_rounds, [] ,epsilon,pdf_prior,noise_prior,xx,0,1);
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
    [xrecon, solutions(:,ILDLC), srecon, iterations(ILDLC),conLDLC]=NBP(H,x',y,sigma, max_rounds, [] ,epsilon,pdf_prior,noise_prior,xx,0,1);
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
%% all went well, we can enter results of this round
for t=1:algs
    results(er,sp,t,1) = results(er,sp,t,1)+times(t);
    results(er,sp,t,3) = results(er,sp,t,3)+iterations(t);
    tic;
    [X_amb L_amb rounds] = Round_and_Local(solutions(:,t),K,H,p,n,y,sigma);
    time = toc;
    avgT = avgT+time;
    totalR = totalR+rounds;
    results(er,sp,t,6) = results(er,sp,t,6)+rounds;
    if x== X_amb(:,1)
        results(er,sp,t,2) = results(er,sp,t,2)+1;
        results(er,sp,t,7) = results(er,sp,t,7)+1;
    end;
    if L(H, p, n, x, y, sigma)>L(H, p, n, X_amb(:,1), y, sigma)
        results(er,sp,t,2) = results(er,sp,t,2)+1;
        results(er,sp,t,4) = results(er,sp,t,4)+1;
    end;
    best(t) = L_amb(1);
end
%% Calc relative rating of algs
for t=1:algs
    results(er,sp,t,5) = results(er,sp,t,5)+sum(unique(sort(best))<=best(t));
    if (sum(unique(sort(best))<=best(t))==1)
        results(er,sp,t,8) = results(er,sp,t,8)+1;
    end
end

end
%% print this rounds (er+sp) results
test = test-skipped;
fprintf('error level     %d     sp level     %d     failed     %d  \n',er,sp,skipped);
fprintf('NBP        : Time- %7.4f     succRate- %6.2f      Iter- %6.2f       MAP- %d       Place- %7.4f       LO- %6.2f\n',results(er,sp,INBP,1)/test, results(er,sp,INBP,2)/test,results(er,sp,INBP,3)/test,results(er,sp,INBP,4),results(er,sp,INBP,5)/test,results(er,sp,INBP,6)/test);
fprintf('IP           : Time- %7.4f     succRate- %6.2f      Iter- %6.2f       MAP- %d       Place- %7.4f       LO- %6.2f\n',results(er,sp,IIP,1)/test, results(er,sp,IIP,2)/test,results(er,sp,IIP,3)/test,results(er,sp,IIP,4),results(er,sp,IIP,5)/test,results(er,sp,IIP,6)/test);
fprintf('CoSaMP: Time- %7.4f     succRate- %6.2f      Iter- %6.2f       MAP- %d       Place- %7.4f       LO- %6.2f\n',results(er,sp,ICoSaMP,1)/test, results(er,sp,ICoSaMP,2)/test,results(er,sp,ICoSaMP,3)/test,results(er,sp,ICoSaMP,4),results(er,sp,ICoSaMP,5)/test,results(er,sp,ICoSaMP,6)/test);
fprintf('GPSR     : Time- %7.4f     succRate- %6.2f      Iter- %6.2f       MAP- %d       Place- %7.4f       LO- %6.2f\n',results(er,sp,IGPSR,1)/test, results(er,sp,IGPSR,2)/test,results(er,sp,IGPSR,3)/test,results(er,sp,IGPSR,4),results(er,sp,IGPSR,5)/test,results(er,sp,IGPSR,6)/test);
fprintf('hardIO: Time- %7.4f     succRate- %6.2f      Iter- %6.2f       MAP- %d       Place- %7.4f       LO- %6.2f\n',results(er,sp,IhIO,1)/test, results(er,sp,IhIO,2)/test,results(er,sp,IhIO,3)/test,results(er,sp,IhIO,4),results(er,sp,IhIO,5)/test,results(er,sp,IhIO,6)/test);
fprintf('NBP-CS : Time- %7.4f     succRate- %6.2f      Iter- %6.2f       MAP- %d       Place- %7.4f       LO- %6.2f\n',results(er,sp,ICS,1)/(test-errors(ICS)), results(er,sp,ICS,2)/(test-errors(ICS)),results(er,sp,ICS,3)/(test-errors(ICS)),results(er,sp,ICS,4),results(er,sp,ICS,5)/(test),results(er,sp,ICS,6)/test);
fprintf('LDLC     : Time- %7.4f     succRate- %6.2f      Iter- %6.2f       MAP- %d       Place- %7.4f       LO- %6.2f\n',results(er,sp,ILDLC,1)/(test-errors(ILDLC)), results(er,sp,ILDLC,2)/(test-errors(ILDLC)),results(er,sp,ILDLC,3)/(test-errors(ILDLC)),results(er,sp,ILDLC,4),results(er,sp,ILDLC,5)/(test),results(er,sp,ILDLC,6)/test);
fprintf('NON        : Time- %7.4f     succRate- %6.2f      Iter- %6.2f       MAP- %d       Place- %7.4f       LO- %6.2f\n',results(er,sp,INON,1)/test, results(er,sp,INON,2)/test,results(er,sp,INON,3)/test,results(er,sp,INON,4),results(er,sp,INON,5)/test,results(er,sp,INON,6)/test);
fprintf('errors: NBP - %d IIP- %d CoSaMP- %d GPSR- %d hardIO- %d NBP-CS- %d\n',errors([1 2 3 4 5 6]));
errors = zeros(1, algs);
if skipped>=10
    assert(false);
end

end
end

if LOFlag
    fprintf('Average time per LO round: %7.4f\n',avgT/totalR);
end
res = results;

%% this is the part of the actual plotting
%%% will work for 1 error level, can choose which, will plot all algorithms
hh=figure;
hold on;
xlabel('Sparce percentage');
ylabel('Success rate');
if LOFlag
    tit = sprintf('Plot of Success rate vs. sparcity, with Local Opt');
    title(tit);
else
    tit = sprintf('Plot of Success rate vs. sparcity, no Local Opt');
    title(tit);
end
for t=1:algs-1
    Xaxis = [0.1:0.1:0.1*Slevels];
    Yaxis = results(er,:,t,2)/test;
    plot(Xaxis,Yaxis,'-x','Color', colors(t), 'LineWidth',2,'DisplayName',names(t,:));
end
hh=plot(Xaxis, results(er,:,algs,2)/test,'--xb', 'LineWidth',2,'DisplayName',names(algs,:));
print -deps /usr0/bickson/fdplot2.eps
saveas(hh,'/usr0/bickson/fdfig2.fig');

%% This is the local opt function
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
save /usr0/bickson/fdres2.mat res;
end
