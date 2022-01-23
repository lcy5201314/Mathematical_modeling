%script for creating figures 3+4 in the paper

% "Fixing the convergence of the GaBP algorithm" by
% J. K. Johnson, D. Bickson and D. Dolev
% In ISIT 2009
% http://arxiv.org/abs/0901.4192
% Code written by Danny Bickson
% May 2009

function[]= run_fix3()
    gamma = [0.055 0.06 0.07 0.08 0.09 0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7]
    %gamma = 0:0.1:5;
%      val = 0.1;
%      gamma = [];
%      for i=1:10
%        val = val*sqrt(2);
%        gamma = [gamma val]; 
%     end
    
    len=length(gamma);
    inner = zeros(1,len);
    outer = inner;
    spec = inner;
    
    for i=1:length(gamma)
        [inner(i),outer(i),spec(i)] = multiuser_detection_fix3(gamma(i));
    end
    figure;
    plot(gamma,inner,'*');
    title('Diagonal weigting vs. inner loop iterations','FontSize', 16);
    ylabel('Avg. number of inner loop iterations','FontSize', 14);
    xlabel('Diagonal weighting \gamma','FontSize', 14);
    figure;
    hold on;
    plot(gamma,outer,'*');
    plot(gamma,outer.*inner,'+');
    title('Diagonal weigting vs. num of iterations','FontSize', 16);
    ylabel('Iterations','FontSize', 14);
    xlabel('Diagonal weighting \gamma','FontSize', 14);
    legend('Outer','Total','FontSize',14);
end