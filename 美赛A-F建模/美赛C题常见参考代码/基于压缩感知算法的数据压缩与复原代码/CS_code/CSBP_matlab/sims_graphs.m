%-------------------
% sims_graphs.m
% Generates graphs for BP paper.
% Dror, 11.25.2008
%-------------------

%-------------------
% plot for L and M varying (sims1)
%-------------------
clear
n=1000;
k=n/10;
SNR=100;
load rep1c_101
E_bp=mean(err_bp_all(:,:,1:ind_r),3); % expected values (can be medians too)
E_m=mean(meas_all(:,:,1:ind_r),3);
xplot_vals=[0.07*n 0.78*n];
subplot(2,1,1)
hold off
plot(E_m(1,1:4),E_bp(1,1:4),'b-*',E_m(2,1:8),E_bp(2,1:8),'b-o',E_m(4,:),E_bp(4,:),'b-x');
legend('L=5','L=10','L=20');
hold on
plot(xplot_vals,ones(1,2)*sqrt(k*SNR),'k--',xplot_vals,ones(1,2)*sqrt(n),'k--');
axis([xplot_vals sqrt(n)*0.5 sqrt(k*SNR)*1.05]);
xlabel('M');
ylabel('MMSE error');
print -depsc error_L_M

%-------------------
% plot for BP vs LP (sims2)
%-------------------
n=500;
k=n/10;
SNR=100;
xplot_vals=[0.07*n 0.78*n];
hold off
load rep2d_101_v5
E_bp=mean(err_bp_all(:,1:ind_r),2); % expected values (can be medians too)
E_lp=mean(err_lp_all(:,1:ind_r),2); 
E_co=mean(err_co_all(:,1:ind_r),2); 
E_m=mean(meas_all(:,1:ind_r),2);
subplot(2,1,1),
plot(E_m,E_co,'b-*',E_m,E_lp,'b-o',E_m,E_bp,'b-x');
legend('CoSaMP','LP','CS-BP');

hold on
plot(xplot_vals,ones(1,2)*sqrt(k*SNR),'k--',xplot_vals,ones(1,2)*sqrt(n),'k--');
axis([xplot_vals sqrt(n)*0.5 sqrt(k*SNR)*1.05]);
xlabel('M');
ylabel('MMSE error');
print -depsc bp_lp_cosamp

%-------------------
% CS-BP with noisy measurements (sims3)
%-------------------
n=1000;
k=n/10;
SNR=100;
load rep3b_80
E_bp=mean(err_bp_all(:,:,1:ind_r),3); % expected values (can be medians too)
E_m=mean(meas_all(:,1:ind_r),2);
xplot_vals=[0.07*n 0.78*n];
subplot(2,1,1)
hold off
plot(E_m,E_bp(:,4),'b-*',E_m,E_bp(:,3),'b-o',E_m,E_bp(:,2),'b-x',E_m,E_bp(:,1),'b-v');
legend('Z=10','Z=5','Z=2','Noiseless');
hold on
plot(xplot_vals,ones(1,2)*sqrt(k*SNR),'k--',xplot_vals,ones(1,2)*sqrt(n),'k--');
axis([xplot_vals sqrt(n)*0.5 sqrt(k*SNR)*1.05]);
xlabel('M');
ylabel('MMSE error');
print -depsc noisy_meas

%-------------------
% CS-BP with wrong model (sims4)
%-------------------
n=1000;
k=n/10;
SNR=100;
load rep4b_103
E_bp=mean(err_bp_all(:,:,1:ind_r),3); % expected values (can be medians too)
E_m=meas_all(:,2,1);
xplot_vals=[0.07*n 0.78*n];
subplot(2,1,1)
hold off
plot(E_m,E_bp(:,5),'b-*',E_m,E_bp(:,4),'b-o',E_m,E_bp(:,3),'b-x',E_m,E_bp(:,2),'b-v');
legend('C=5','C=4','C=3','C=2');
hold on
plot(xplot_vals,ones(1,2)*sqrt(k*SNR),'k--',xplot_vals,ones(1,2)*sqrt(n),'k--');
axis([xplot_vals sqrt(n)*0.5 sqrt(k*SNR)*1.05]);
xlabel('M');
ylabel('MMSE error');
print -depsc wrong_mixture
