function [xrecon, mrecon, srecon, pdf, pdf_xx, pdf_prior]=decoder(x,measvec,n,k,l,r,sig_1,sig_0,sigma_Z,iter,phi,phisign,aux,model_order,self_indexN,self_indexM,dispind,boundx,epsilon,gamma_pdbp,gamma_mdbpf,gamma_mdbpb);

%---------------
% initializations
%---------------
b1=1;
b2=1;
b3=1;
n=length(x);      
m=length(measvec);
xrecon=zeros(1,n); % max likelihood
mrecon=zeros(1,n); % average 
srecon=zeros(1,n); % std_dev
pdf_N_to_M=zeros(n, r, model_order);
pdf_M_to_N=zeros(m, l, model_order);
xx=1:model_order;  xx=xx-(model_order+1)/2;  xx=xx/max(xx);
xx=xx*boundx; % values over which pdf is sampled
pdf_xx=xx;
delta=xx(2)-xx(1);

%---------------
% compute signal node prior
%---------------
pdf_prior= (k/n)*normpdf(xx,0,sig_1);
if (sig_0 > epsilon)
  pdf_prior= pdf_prior  +    (1-k/n)*normpdf(xx,0,sig_0);
else
  in2=find(abs(xx)<epsilon);
  pdf_prior(in2)=pdf_prior(in2) + (1-k/n);
end
pdf_prior=pdf_prior/sum(pdf_prior);

y_noise=normpdf(xx,0,sigma_Z);	% noise prior
y_noise=y_noise/sum(y_noise);

pdf_prior=ifftshift(pdf_prior);
pdf_xx=ifftshift(pdf_xx);
xx=ifftshift(xx);
y_noise=ifftshift(y_noise);
pdf=zeros(n,length(xx));

%---------------
% BP ITERATIONS 
%---------------
for it=1:iter	
	%---------------
   % FORWARD ITERATION - from measurements to signal
	%---------------
   for i=1:n
      if (it==1) % initial pdf
         for rr=1:r
            pdf_N_to_M(i, rr, :)=pdf_prior(:);
         end
      else         
         neighbors=aux(:,i)';
         neighbors=setdiff_shri(neighbors,0);
         ln=length(neighbors);
         self_index=self_indexN(i,:);
         %---------------
         % ESTIMATE SIGNAL COEFF x(i)
         %---------------
         pdf_res=[];
         for jj=1:ln 
            pdf_tmp=reshape(pdf_M_to_N(neighbors(jj),self_index(jj),:),1, model_order);
            pdf_tmp=pdf_tmp+epsilon;
            if (length(pdf_res)==0) % first time
               pdf_res=(pdf_tmp);
            else % convolution step
               pdf_res=mulpdf(pdf_res,pdf_tmp);
            end
         end
         pdf_res=mulpdf(pdf_res, pdf_prior);
         if (it>1) % damping
            pdf_res=gmean(pdf_res,pdf(i,:),gamma_pdbp,0); %gamma_pdbp=0: true BP
         end
         [mtt, stt, maxtt]=meanvarmaxpdf((pdf_res), xx); % computes statistics
         xrecon(i)=maxtt; % max likelihood
         mrecon(i)=mtt; % mean
         srecon(i)=stt; % std_dev around mean
         pdf(i,:)=pdf_res; % store result
         % END OF ESTIMATE OF SIGNAL
         for j=1:ln   % to send next message
            pdf_tmp=reshape(pdf_M_to_N(neighbors(j),self_index(j),:),1, model_order);
            pdf_tmp=pdf_tmp+epsilon; % for stability
            [pdf_tosend]=divpdf(pdf_res,pdf_tmp);
            pdf_tosend=pdf_tosend/sum(pdf_tosend);
            tmptmp=pdf_tosend;
            if (it>1) % MDBPF
               prevpdf=reshape(pdf_N_to_M(i, j,:),1,model_order);;
               tmptmp=gmean(pdf_tosend,prevpdf,gamma_mdbpf,0);
            end
            pdf_N_to_M(i, j, :)=tmptmp/sum(tmptmp);
         end %of for j
      end
   end
   if (it >=2) % display and break on last iteration
      dispvec_anderrors(mrecon, phi, phisign, measvec, dispind, x);
      if (it==iter)
         break;
      end
   end
	%---------------
   % BACKWARD ITERATION - from signal to measurements
	%---------------
   for i=1:m
      neighbors=phi(i,:);
      neighbors=setdiff_shri(neighbors,0);
      ln=length(neighbors);
      phisigni=phisign(i,1:ln);
      self_index=self_indexM(i,:);
      pdf_res_all=[];
      for jj=1:ln % process neighbors
         tmptmp=[reshape(pdf_N_to_M(neighbors(jj),self_index(jj),:),1, model_order)];
         if ((phisigni(jj) < 0) & (b1))
            tmptmp=reverse(tmptmp,model_order);
         end
         if (length(pdf_res_all)==0) % first time
            pdf_res_all=(fft((tmptmp)));
            pdf_res_all=pdf_res_all+epsilon;
         else           
            pdf_res_all=convpdf_fft(pdf_res_all, tmptmp, epsilon);
         end
      end
      if (sigma_Z>epsilon)
         [pdf_res_all]=convpdf_fft(pdf_res_all, y_noise, epsilon);
      end
      for j=1:ln   %To send next message
         tmptmp=[reshape(pdf_N_to_M(neighbors(j),self_index(j),:),1, model_order)];
         mvi=measvec(i);
         if ((phisigni(j) < 0) & (b2))
            tmptmp=reverse(tmptmp,model_order);
         end
         pdf_res=unconvpdf_fft(pdf_res_all, tmptmp, epsilon);
         pdf_res=abs((ifft((pdf_res))));
         pdf_res=shiftpdf_fft(pdf_res, mvi, delta, model_order);  
         pdf_res=pdf_res/sum(pdf_res);
         if ((phisigni(j) < 0) & (b3)) % MDPB
            pdf_res=reverse(pdf_res,model_order);
         end
         tmptmp=pdf_res;
         if (it>1) % MDBPF
            prevpdf=reshape(pdf_M_to_N(i, j, :), 1, model_order);;
            tmptmp=gmean(pdf_res, prevpdf, gamma_mdbpb,0);
         end
         pdf_M_to_N(i, j, :)=tmptmp/sum(tmptmp);
      end
   end
end
