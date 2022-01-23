%Code written by Harel Avissar and Danny Bickson
%Suplamentary material for the paper: 
%Distributed fault identification via non-parametric belief propagation.
%By D. Bickson, H. Avissar, D. Dolev, S. P. Boyd, A. Ihler and D. Baron.
%Submitted for publication. Aug 2009.
%Code available from http://www.cs.huji.ac.il/labs/danss/p2p/gabp/

%code for testing SDP formulation of the fault det problem.
%requires CVX to run. http://www.stanford.edu/~boyd/cvx/

clear all

avgT = 0;
totalR = 0;
LOFlag = 0;
K = 10;
n = 100;%15;
m = 50;%10;
Slevels = 4; %sparcity levels to check
Elevels = 1; %error levels to check
Sstart = 1;
Estart = 1;
algs = 9; %number of algorithems compare
params = 8; %number of parameters to monitor on each alg
% params are : 1- time      2- success rate         3- iterations to
% converge      4- number of times when found MAP solution better than real
% solution      5- ranking of solution          6- number of LO rounds
% 7- number of exact solution           8- number of times on top
INBP = 1; IIP = 2; ICoSaMP = 3; IGPSR = 4; IhIO = 5; ICS = 6; INON = 7; ILDLC = 8; ISDP = 9;%index for each alg
names = ['NBP   ';'IP    ';'CoSaMP';'GPSR  ';'hardIO';'CS    ';'NON   ';'LDLC  ';'SDP   '];
colors = ['b', 'r', 'c', 'm', 'k', 'g', 'y', 'w', 'w'];
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
test =125;
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
    
    
%% cvx for SDP
tic;
cvx_quiet(true);
lambda = log((1-p)/p)*ones(n,1);
 cvx_begin
     variable g(n)
     minimize((1/(2*sigma^2))*square_pos(norm(A*g-Y,2))+lambda'*g)
     g >= 0
     g <= 1
 cvx_end
 solutions(:,ISDP) = 2*(g-0.5);
 iterations(ISDP) = 0;
 times(ISDP) = toc;
 %l_min = cvx_optval-(1/(2*sigma^2))*norm(y)^2
 %cvx_optval
 %x
    
%% all went well, we can enter results of this round
for t=algs:algs
    results(er,sp,t,1) = results(er,sp,t,1)+times(t);
    results(er,sp,t,3) = results(er,sp,t,3)+iterations(t);
    tic;
    [X_amb L_amb rounds] = Round_and_Local(solutions(:,t),K,H,p,n,y,sigma,LOFlag);
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
% for t=1:algs
%     results(er,sp,t,5) = results(er,sp,t,5)+sum(unique(sort(best))<=best(t));
%     if (sum(unique(sort(best))<=best(t))==1)
%         results(er,sp,t,8) = results(er,sp,t,8)+1;
%     end
% end

%% Load old results for the others
load('results.mat');
if LOFlag
    results(:,:,1:algs-1,:) = resLO;
else
    results(:,:,1:algs-1,:) = res;
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
fprintf('CVX        : Time- %7.4f     succRate- %6.2f      Iter- %6.2f       MAP- %d       Place- %7.4f       LO- %6.2f\n',results(er,sp,ISDP,1)/test, results(er,sp,ISDP,2)/test,results(er,sp,ISDP,3)/test,results(er,sp,ISDP,4),results(er,sp,ISDP,5)/test,results(er,sp,ISDP,6)/test);
fprintf('errors: NBP - %d IIP- %d CoSaMP- %d GPSR- %d hardIO- %d NBP-CS- %d\n',errors([1 2 3 4 5 6]));
fprintf('NBP        : Excact- %d \n',results(er,sp,INBP,7));
fprintf('IP           : Excact- %d \n',results(er,sp,IIP,7));
fprintf('CVX        : Excact- %d \n',results(er,sp,ISDP,7));
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
for t=1:algs-2
    Xaxis = [0.1:0.1:0.1*Slevels]*100;
    Yaxis = (results(er,:,t,2)/test)*100;
    plot(Xaxis,Yaxis,'-x','Color', colors(t), 'LineWidth',2,'DisplayName',names(t,:));
end
plot(Xaxis, (results(er,:,algs-1,2)/test)*100,'--xb', 'LineWidth',2,'DisplayName',names(algs-1,:));
plot(Xaxis, (results(er,:,algs,2)/test)*100,'--xr', 'LineWidth',2,'DisplayName',names(algs,:));
