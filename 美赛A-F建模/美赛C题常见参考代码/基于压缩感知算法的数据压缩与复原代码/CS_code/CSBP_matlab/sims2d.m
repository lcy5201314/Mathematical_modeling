%-------------------
% sims2.m
% Runs second set of sims for BP paper.
% For different M and fixed L, compare LP and BP.
% Dror, 11.25.2008
%-------------------
% version d: also uses cosamp, Dror, 12.9.2008
%-------------------

clear
rand('state',sum(100*clock));

%-------------------
% Parameter definitions
%-------------------
restart=1;	% 1 means start sims from scratch else read in old file
%Signal and noise parameters
num_reps=101;
iter=15; % Number of iterations in BP
n=500; %Signal length
k=n/10; %Sparsity
l=20; % # nonzeros per row

%%NOTE: *actual* number of msmts is the no. rows in phi below
%m_all=round(n*(0.1:0.1:0.8)); 
r_all=1:17;

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
len_r=length(r_all);
if (restart==1) % initialize values
   ind_r=0;
   err_bp_all=zeros(len_r,num_reps);
   err_lp_all=zeros(len_r,num_reps);
   err_co_all=zeros(len_r,num_reps);
   meas_all=zeros(len_r,num_reps);
else % load previously stored values
   load rep1d_5
end;
for ind_r=ind_r+1:num_reps
   fprintf('BIG ITERATION NUMER %4d\n',ind_r);
   for ind_R=1:len_r
      R=r_all(ind_R);
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
         err_bp_all(ind_R,ind_r)=norm(mrecon-x);
         meas_all(ind_R,ind_r)=num_rows;
         fprintf('Time=%8.2f\n\n',cputime-t1);
         %-------------------
         % LP
         %-------------------
         % initial guess = min energy
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
         err_lp_all(ind_R,ind_r)=norm(xh'-x);
         fprintf('Time=%8.2f\n\n',cputime-t1);
         %-------------------
         % cosamp
         %-------------------
         fprintf('Solving cosamp\n');
         xco=cosamp(measvec,phi2,k,100);
         fprintf('Cosamp error: %6.2f\n',norm(xco(:,100)-x'));
         err_co_all(ind_R,ind_r)=norm(xco(:,100)-x');
         fprintf('Time=%8.2f\n\n',cputime-t1);
      end;
   end;
   %-------------------
   % save file and graphics
   %-------------------
   eval(['save rep2d_' num2str(ind_r) ' ind_r meas_all err_bp_all err_lp_all err_co_all']);
   E_bp=mean(err_bp_all(:,1:ind_r),2); % expected values (can be medians too)
   E_lp=mean(err_lp_all(:,1:ind_r),2); 
   E_co=mean(err_co_all(:,1:ind_r),2); 
   E_m=mean(meas_all(:,1:ind_r),2);
   xplot_vals=[0.07*n 0.83*n];
   hold off
   plot(E_m,E_bp,'b-*',E_m,E_lp,'b-o',E_m,E_co,'b-v');
   legend('BP','LP','CoSaMP');
   hold on
   plot(xplot_vals,ones(1,2)*sqrt(k*SNR),'g--',xplot_vals,ones(1,2)*sqrt(n),'g--');
   axis([xplot_vals sqrt(n)*0.5 sqrt(k*SNR)*1.05]);
   xlabel('M');
   ylabel('Reconstruction error');
   eval(['print -depsc sim2d_' num2str(ind_r)]);
end;
