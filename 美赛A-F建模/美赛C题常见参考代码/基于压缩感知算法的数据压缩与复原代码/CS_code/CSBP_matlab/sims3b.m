%-------------------
% sims3a.m
% For fixed M and fixed L, run CS-BP on different noise levels.
% Dror, 12.9.2008
%-------------------
% sims3b: modified noise levels for better graphics.
%-------------------

clear
rand('state',sum(100*clock));

%-------------------
% Parameter definitions
%-------------------
restart=0;	% 1 means start sims from scratch else read in old file
%Signal and noise parameters
num_reps=301;
iter=15; % Number of iterations in BP
n=1000; %Signal length
k=n/10; %Sparsity
l=20; % # nonzeros per row

%%NOTE: *actual* number of msmts is the no. rows in phi below
%m_all=round(n*(0.1:0.1:0.8)); 
r_all=1:17;

SNR=100;
sigma_signal=sqrt(SNR); 
sigma_noise=1;  %Noise IN the signal x

%Damping for Belief Prop
gamma_mdbpf=0.35;
gamma_mdbpb=0.35;
gamma_pdbp=0.0;

model_order=243	% Number of sampling points 
% IMORTANT: FFT runs faster when this value can be partitioned nicely

%-------------------
% Noise levels
% This is the novelty in this file
%-------------------
noise_all=[1e-20 2 5 10];
%s_n=1e-20;		%Noise in the measurements y (noisy measurements)

%-------------------
% MAIN LOOP
%-------------------
t1=cputime;
len_r=length(r_all);
len_noise=length(noise_all);
if (restart==1) % initialize values
   ind_r=0;
   err_bp_all=zeros(len_r,len_noise,num_reps);
   meas_all=zeros(len_r,len_noise,num_reps);
else % load previously stored values
   load rep3b_72
%   load rep3b_43
%   ind_r2=ind_r;
%   err_bp2=err_bp_all;
%   meas2=meas_all;
%   load rep3b_16
%   err_bp_all(:,:,ind_r+1:ind_r+ind_r2)=err_bp2(:,:,1:ind_r2);
%   meas_all(:,:,ind_r+1:ind_r+ind_r2)=meas2(:,:,1:ind_r2);
%   ind_r=ind_r+ind_r2, pause
%   clear ind_r2 meas2 err_bp2
end;
for ind_r=ind_r+1:num_reps
   fprintf('BIG ITERATION NUMER %4d\n',ind_r);
   for ind_R=1:len_r
      R=r_all(ind_R);
      if (R)
         %-------------------
         % ITERATE OVER NOISE LEVELS
         %-------------------
         for ind_noise=1:len_noise
%             s_n=noise_all(ind_noise);
				ynoise_sigma_actual=noise_all(ind_noise); %used in driver.m; the noise in measurments
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
            err_bp_all(ind_R,ind_noise,ind_r)=norm(mrecon-x);
            meas_all(ind_R,ind_noise,ind_r)=num_rows;
            fprintf('Time=%8.2f\n\n',cputime-t1);
         end;
      end;
   end;
   %-------------------
   % save file and graphics
   %-------------------
   eval(['save rep3b_' num2str(ind_r) ' ind_r meas_all err_bp_all']);
   E_bp=mean(err_bp_all(:,:,1:ind_r),3); % expected values (can be medians too)
   E_m=mean(meas_all(:,1:ind_r),2);
   xplot_vals=[0.07*n 0.83*n];
   hold off
   plot(E_m,E_bp(:,4),'b-*',E_m,E_bp(:,3),'b-o',E_m,E_bp(:,2),'b-v',E_m,E_bp(:,1),'b-p');
   legend('N=10','N=5','N=2','Noiseless');
   hold on
   plot(xplot_vals,ones(1,2)*sqrt(k*SNR),'k--',xplot_vals,ones(1,2)*sqrt(n),'k--');
   axis([xplot_vals sqrt(n)*0.5 sqrt(k*SNR)*1.05]);
   xlabel('M');
   ylabel('MMSE error');
   eval(['print -depsc sim3b_' num2str(ind_r)]);
end;
