%-------------------
% sims1.m
% Runs first set of sims for BP paper.
% Different M and L for BP.
% Based on Shri's original main.m script.
% Dror, 11.21.2008
%-------------------

clear
rand('state',sum(100*clock));

%-------------------
% Parameter definitions
%-------------------
restart=0;	% 1 means start sims from scratch else read in old file
%Signal and noise parameters
num_reps=101;
iter=15; % Number of iterations in BP
n=1000; %Signal length
k=n/10; %Sparsity

% Parameters simulated with different values
l_all=[5 10 15 20];

%%NOTE: *actual* number of msmts is the no. rows in phi below
%m_all=round(n*(0.1:0.1:0.8)); 
r_all=[1:4 zeros(1,10);1:8 zeros(1,6);2:11 zeros(1,4); 1:14]

SNR=100;
sigma_signal=sqrt(SNR); 
sigma_noise=1;  %Noise IN the signal x
s_n=1e-20;		%Noise in the measurements y (noisy measurements)
ynoise_sigma_actual=s_n; %used in driver.m; the noise in measurments

%Damping for Belief Prop
gamma_mdbpf=0.35;
gamma_mdbpb=0.35;
gamma_pdbp=0.0;

model_order=243	% Number of sampling points 
% IMORTANT: FFT runs faster when this value can be partitioned nicely

%-------------------
% MAIN LOOP
%-------------------
t1=cputime;
len_l=length(l_all);
len_r=length(r_all);
if (restart==1) % initialize values
   ind_r=0;
   err_bp_all=zeros(len_l,len_r,num_reps);
   err_lp_all=zeros(len_l,len_r,num_reps);
   meas_all=zeros(len_l,len_r,num_reps);
else % load previously stored values
   load rep1c_5
end;
for ind_r=ind_r+1:num_reps
   fprintf('BIG ITERATION NUMER %4d\n',ind_r);
   for ind_l=1:len_l
      for ind_R=1:len_r
         l=l_all(ind_l);
         R=r_all(ind_l,ind_R);
         if (R)
            %-------------------
            % Generate signal
            %-------------------
            disp ('GENERATING THE SIGNAL X...');
            [x, heavyind]=generatex_noisy(n, k, sigma_signal, sigma_noise);
            disp(sprintf('l2 norm of x: %g', norm(x) )); 
            %-------------------
            % Generate measurement matrix
            %-------------------
            [phi]=gen_phi(n, l,R);
            phisign=randn(size(phi));  phisign=sign(phisign);
            num_rows=size(phi,1);
            phi2=zeros(num_rows,n);
            for ind_row=1:num_rows
               phi2(ind_row,phi(ind_row,:))=phisign(ind_row,:);
            end;
            %-------------------
            % BP decoding
            %-------------------
            driver
            err_bp_all(ind_l,ind_R,ind_r)=norm(mrecon-x);
            meas_all(ind_l,ind_R,ind_r)=num_rows;
            %-------------------
            % LP
            %-------------------
            % initial guess = min energy
            if (0)
               x0=phi2'*inv(phi2*phi2')*measvec;
               x0p = [max(x0,0); -min(x0,0)];
               c = ones(2*n,1);
               lower = zeros(2*n,1);
               upper = 100*ones(2*n,1);
               A_lp = [phi2 -phi2];
               % solve the LP
               fprintf('Solving LP\n');
               xhlp = lp_pdco(c, measvec, A_lp, x0p, lower, upper, 150000, 1e-4);
               Ntmp=length(xhlp) / 2;
               xh = xhlp(1:Ntmp) - xhlp(Ntmp+1:2*Ntmp);
               fprintf('LP error: %6.2f\n',norm(xh'-x));
               err_lp_all(ind_l,ind_R,ind_r)=norm(xh'-x);
            end;
            fprintf('Time=%8.2f\n\n',cputime-t1);
         end;
      end;
      eval(['save rep1c_temp ind_r meas_all err_bp_all err_lp_all']);
   end;
   %-------------------
   % save file and graphics
   %-------------------
   eval(['save rep1c_' num2str(ind_r) ' ind_r meas_all err_bp_all err_lp_all']);
   E_bp=mean(err_bp_all(:,:,1:ind_r),3); % expected values (can be medians too)
   E_m=mean(meas_all(:,:,1:ind_r),3);
   xplot_vals=[0.07*n 0.83*n];
   hold off
   plot(E_m(1,1:4),E_bp(1,1:4),'b-*',E_m(2,1:8),E_bp(2,1:8),'b-o',...
      E_m(3,1:10),E_bp(3,1:10),'b-x',E_m(4,:),E_bp(4,:),'b-v');
   legend('L=5','L=10','L=15','L=20');
   hold on
   plot(xplot_vals,ones(1,2)*sqrt(k*SNR),'g--',xplot_vals,ones(1,2)*sqrt(n),'g--');
   axis([xplot_vals sqrt(n)*0.5 sqrt(k*SNR)*1.05]);
   xlabel('M');
   ylabel('Reconstruction error');
   eval(['print -depsc sim1c_' num2str(ind_r)]);
end;
