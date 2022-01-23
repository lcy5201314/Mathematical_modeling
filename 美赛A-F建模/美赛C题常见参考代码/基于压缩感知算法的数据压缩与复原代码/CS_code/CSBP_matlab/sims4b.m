%-------------------
% sims4a.m
% For fixed L, run CS-BP on mixture Gaussian and compressble sources.
% Dror, 12.9.2008
%-------------------
% sims4b.m: instead of compressible, comparing 2-, 3-, 4-state mixtures.
% Dror, 12.11.2008
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
S=0.1; % sparsity rate
k=n*S; % # big 
l=20; % # nonzeros per row

%%NOTE: *actual* number of msmts is the no. rows in phi below
%m_all=round(n*(0.1:0.1:0.8)); 
r_all=1:17;

SNR=100;
sigma_signal=sqrt(SNR); 
sigma_noise=1;  %Noise IN the signal x
ynoise_sigma_actual=1e-20; %used in driver.m; the noise in measurments

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
   err_bp_all=zeros(len_r,2,num_reps);
   meas_all=zeros(len_r,2,num_reps);
else % load previously stored values
   load rep4b_41
   ind_r2=ind_r;
   err_bp2=err_bp_all;
   meas2=meas_all;
   load rep4b_54
   err_bp_all(:,:,ind_r+1:ind_r+ind_r2)=err_bp2(:,:,1:ind_r2);
   meas_all(:,:,ind_r+1:ind_r+ind_r2)=meas2(:,:,1:ind_r2);
   ind_r=ind_r+ind_r2, pause
   clear ind_r2 meas2 err_bp2
end;
for ind_r=ind_r+1:num_reps
   fprintf('BIG ITERATION NUMER %4d\n',ind_r);
   for ind_R=1:len_r
      R=r_all(ind_R);
      if (R)
         %-------------------
         % TWO DIFFERENT SIGNALS
         %-------------------
         for mixture_num=2:5
            %-------------------
            % Generate signal
            %-------------------
            fprintf('Number of Mixture components %2d\n',mixture_num);
            sigma_signal=sqrt(SNR-sigma_noise^2); 
            relative_amplitude=1:mixture_num-1; % generate amplitudes
            relative_amplitude=relative_amplitude/norm(relative_amplitude);
            disp ('GENERATING THE SIGNAL X...');
            % [x, heavyind]=generatex_noisy(n, k, sigma_signal, sigma_noise);
            x=randn(1,n)*sigma_noise; % generate small coeffs
            loc1=ceil(rand(1,n)/S); % generate location of bigs
            for ind2=mixture_num-1:-1:1 % iterate over component num
               loc_here=find(loc1==ind2);
               x(loc_here)=x(loc_here)+randn(1,length(loc_here))*sigma_signal*relative_amplitude(ind2);
            end;
            [tmp,heavyind]=sort(abs(x));
            heavyind=heavyind(n-9:n);
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
            err_bp_all(ind_R,mixture_num,ind_r)=norm(mrecon-x);
            meas_all(ind_R,mixture_num,ind_r)=num_rows;
            fprintf('Time=%8.2f\n\n',cputime-t1);
         end;
      end;
   end;
   %-------------------
   % save file and graphics
   %-------------------
   eval(['save rep4b_' num2str(ind_r) ' ind_r meas_all err_bp_all']);
   E_bp=mean(err_bp_all(:,:,1:ind_r),3); % expected values (can be medians too)
   E_m=meas_all(:,2,1);
   xplot_vals=[0.07*n 0.83*n];
   hold off
   plot(E_m,E_bp(:,5),'b-*',E_m,E_bp(:,4),'b-o',E_m,E_bp(:,3),'b-h',E_m,E_bp(:,2),'b-v');
   legend('C=5','C=4','C=3','C=2');
   hold on
   plot(xplot_vals,ones(1,2)*sqrt(k*SNR),'g--',xplot_vals,ones(1,2)*sqrt(n),'g--');
   axis([xplot_vals sqrt(n)*0.5 sqrt(k*SNR)*1.05]);
   xlabel('M');
   ylabel('MMSE error');
   eval(['print -depsc sim4b_' num2str(ind_r)]);
end;
