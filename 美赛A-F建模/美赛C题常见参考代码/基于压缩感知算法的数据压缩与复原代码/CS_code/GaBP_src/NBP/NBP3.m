%-----------------------------------------------
% Original code by Shriram Sarvotham, 2006. 
% Cleaned up a bit by Dror Baron, November 2008.
% Full original code is found on: http://www.ece.rice.edu/~drorb/CSBP/
%-----------------------------------------------
% Code rewritten by Danny Bickson, March 2009.
% Impvoed accuracy of computation, added support for non binary matrices 
% Added support for arbitrary self potentials
% New code is avaialble from
% http://www.cs.huji.ac.il/labs/danss/p2p/gabp/index.html
%------------------------------------------------
% Nonparametric belief propagation implementation
% Channel model is: y = Ax + w, where w is AWGN noise ~ N(0,sigma_Z)
% priors of x and y are Gaussian mixtures
% Input:
% A - m x n - real linear transformation matrix
% x - n x 1 - hidden vector
% y - m x 1 - observation vector
% sigma_Z - noise level
% max_iter - maximal number of iterations
% dispind - a vector indices to display
% epsilon - small value to add to fft to avoid division by zero
% pdf_prior - prior of x
% noise_prior - prior of y
% xticks - points to sample the Gaussian mixture


function [xrecon, mrecon, srecon]=NBP(A,x,y,sigma_Z,max_iter,...
    dispind,epsilon,pdf_prior,noise_prior,xticks);

n=length(x);      
m=length(y);
M2N = cell(m,n);
N2M = cell(m,n);
model_order = length(xticks);
delta=xticks(2)-xticks(1);
pdf_prior=ifftshift(pdf_prior);
xticks=ifftshift(xticks);
noise_prior=ifftshift(noise_prior);
last_norm = inf;

%---------------
% BP ITERATIONS 
%---------------
for it=1:max_iter	

   %---------------
   % FORWARD ITERATION - from signal to measurement
   % Calculate the product of all incoming mixtures except the one the
   % message is sent to.
   % The product is simply point by point product
   %---------------
   for i=1:n
       
       % For each neighbor of i
      neighbors = find(A(:,i)~=0)';
      ln=length(neighbors); 
      
      if (it==1) % initial round - send the signal prior
         for j=1:ln
            N2M{neighbors(j),i} = pdf_prior(:)';
         end
      else  % round >= 2        
         for j=1:ln 
            m_ji = M2N{neighbors(j),i} + epsilon;
            verify_pdf(m_ji);
            if (j == 1) % first time
               pdf_res=m_ji;
            else 
               pdf_res=mulpdf(pdf_res,m_ji);
            end
         end
         pdf_res=mulpdf(pdf_res, pdf_prior);
         [mrecon(i), srecon(i),xrecon(i)]=meanvarmaxpdf((pdf_res), xticks); % computes statistics
         for j=1:ln   % to send next message
            m_ji = M2N{neighbors(j),i}+epsilon;
            N2M{neighbors(j),i}= divpdf(pdf_res,m_ji);
         end
      end
   end
   if (it >=2) % display and break on last iteration
      [flag,last_norm] = dispvec_anderrors_gabp(mrecon, A, y, dispind, x, last_norm);
      if (it==max_iter)
         break;
      end
   end
   %---------------
   % BACKWARD ITERATION - from measurement to signal
   % Calculate convolution - which is a product in the FFT domain
   %---------------
   for i=1:m
      neighbors=find(A(i,:)~=0);
      ln=length(neighbors);
      for j=1:ln % process neighbors
         m_ij = N2M{i,neighbors(j)};
         m_ij = pdf_integral(m_ij, model_order, A(i,neighbors(j)),xticks,delta,epsilon);
         if (j == 1) % first time
            %pdf_res_all=(fft([m_ij zeros(1,model_order - 1)]))+epsilon;
            pdf_res_all=m_ij;
            pdf_res_all2=(fft([m_ij]))+epsilon;
         else           
            %pdf_res_all=mulpdf_fft(pdf_res_all, [m_ij zeros(1,model_order - 1)], epsilon);
            pdf_res_all=conv(pdf_res_all,m_ij);
            pdf_res_all2=mulpdf_fft(pdf_res_all2, [m_ij], epsilon);
         end
      end
      
      %convolve with the self potential of the noise
      if (sigma_Z>epsilon)
         %pdf_res_all=mulpdf_fft(pdf_res_all, [noise_prior zeros(1,model_order - 1)], epsilon);
         pdf_res_all=conv(pdf_res_all,noise_prior);
         pdf_res_all2=mulpdf_fft(pdf_res_all2, [noise_prior], epsilon);
      end
      
      for j=1:ln   %To send next message
         m_ij = N2M{i,neighbors(j)};
         m_ij = pdf_integral(m_ij, model_order, A(i,neighbors(j)),xticks,delta,epsilon);
         %unconvolve with the message node from j to i
         %pdf_res=divpdf_fft(pdf_res_all, [m_ij zeros(1,model_order - 1)], epsilon);
         pdf_res=deconv(pdf_res_all, m_ij);
         pdf_res2=divpdf_fft(pdf_res_all2, [m_ij], epsilon);
         % get back from the FFT to the real domain
         %pdf_res=abs((ifft((pdf_res))));
         %pdf_res=ifftshift(pdf_res);
         pdf_res2=abs((ifft((pdf_res2))));
         % compute y(i) - current mixture
         %pdf_res=shiftpdf_fft(pdf_res, y(i), delta, model_order);  
         %l=(model_order-1)/2;
         %pdf_res = pdf_res(l:l+model_order-1);
         %pdf_res = fftshift(pdf_res);
         pdf_res=shiftpdf_fft_gabp(pdf_res, y(i), delta, xticks);  
         %pdf_res=fftshift(pdf_res);
         pdf_res2=shiftpdf_fft_gabp(pdf_res2, y(i), delta, xticks);  
         M2N{i,neighbors(j)} = pdf_integral(pdf_res, model_order, A(i,neighbors(j)),xticks,delta,epsilon);
         verify_pdf(M2N{i,neighbors(j)});
     end
   end
   
 
   
end
% Handle computation of integral.
% For edges with weight 1 - does not do anything
% For edges with weight -1 - reverses the mixture
% For real non zero edges, calculates interpolation
function [npdf]= pdf_integral(pdf, model_order, edge,xticks,delta,epsilon)
   assert(edge ~= 0);
    if (edge == -1)
       npdf=reverse_gabp(pdf,model_order);
   elseif  (edge == 1)
       npdf = pdf;
    else
       %old_mean = fft_max(xticks,pdf);
       old_mean = sum(xticks.*pdf);
       npdf= interp1(xticks, pdf, xticks./abs(edge));
       npdf(isnan(npdf)) = epsilon; 
       if (edge < 0)
           npdf = reverse_gabp(npdf, model_order);
       end
           
       npdf=npdf./sum(npdf);
       new_mean = sum(xticks.*npdf);
       if (length(new_mean) == 1 && (abs(new_mean - old_mean * edge) >= 2*delta))
           %assert(abs(new_mean - old_mean * edge) < delta);%DB TODO 
       end
   end
end
% point by point multiplication
function c=mulpdf(a,b);
c=a.*b;
c=c/sum(c);
end
% point by point division
function c=divpdf(a,b);
c=a./b;
c=c/sum(c);
end
% point by point multiplication in the FFT domain
 function c=mulpdf_fft(a,b,epsilon)
tmp=fft((b));
tmp=max(epsilon,tmp);
c=mulpdf(a,tmp);
 end
 % point by point division in the FFT domain
function c=divpdf_fft(a,b,epsilon)
tmp=fft((b));
tmp=max(epsilon,tmp);
c=divpdf(a,tmp);
end
% sanity checks
function []=verify_pdf(pdf);
   if (sum(pdf)<0 || sum(pdf) > 1.1)
    assert(sum(pdf) > 0 && sum(pdf)<=1);
   end
   assert(sum(isnan(pdf)) == 0);
   assert(sum(~isreal(pdf)) == 0);
end

function [m,sig,mp]=meanvarmaxpdf(pdf, xx);
    pdf=pdf/sum(pdf);
    [mm,ind]=max(pdf);
    mp=xx(ind);
    m=sum(xx.*pdf);
    sig=sqrt(sum(xx.*xx.*pdf)-m*m);
end

function [flag,tnorm] =dispvec_anderrors_gabp (v, A, y, dispind, x, last_norm);
flag = 0;
l=length(dispind);
s='[';
for i=1:l
  s=sprintf('%s %7.2f',s,v(dispind(i)));
end
s=sprintf('%s]',s);
mv=A*v';
er=norm(y-mv);
tnorm = norm(x-v);
s=sprintf('%s (%7.4f) (%7.4f)',s, er, tnorm);
disp(s);
if (last_norm < tnorm)
    flag  = 1;
end
end

% DB: Since the array is fftshifted, the first position is zero (odd sample size assumed), no need to
% reverse it. The other values are reversed, which means we computed
% integral with edge weight of -1.
function op=reverse_gabp(ip, m)
    op=[ip(1) ip(m:-1:2)];
end


%move the new mean to "Offset - old mean"
function pdf=shiftpdf_fft_gabp(pdf, offset, delta, xx)

    %old_pos = fft_max(xx,pdf);
    old_pos = sum(xx.*pdf);
    pdf = fliplr(pdf);
    pdf = interp1(xx,pdf,xx-offset);
 
    pdf(isnan(pdf)) = 0;
    pdf = pdf./sum(pdf);
    %new_pos = fft_max(xx,pdf);
    new_pos = sum(xx.*pdf);
    if (length(old_pos) == 1 && (abs((offset - old_pos) - new_pos) >= 6*delta))
           %assert(abs((offset - old_pos) - new_pos) < delta); %TODO DANNY
    end
    %pdf = fftshift(pdf);
end

end