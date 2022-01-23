
% Non-parametric belief propagation code
% Written by Danny Bickson.
% updated: 3-March-2008
%
%
% input: H - inverse covariance sqaure matrix size n x n containing edge weights. H_ij is inverse 
% covariance between nodes i and j
% self_pot - an array of kde() objects containing the self potentials. The
% self potentials are Gaussian mixtures
% max_iter - maximal number of iterations
% epsilon - convergence threshold
% output: x - a kde of size nx1 containing the beliefs
%

function [meanx,varx,maxx,xx]=UNBP(H,x,y,self_pot,max_iter,dispind, epsilon,xticks)
   
    %options = [-1 1]; % the allowed integers
    %density = 30; % max number of Gaussians in each mixture
    
    model_order = length(xticks);
    delta=xticks(2)-xticks(1);

    n = size(H,1); % number of check nodes
    assert(size(H,2) == size(H,1));
    
    % variable messages
    M = cell(n,n);
    old_M = cell(n,n);
    M_shift = zeros(n,n);
    old_M_shift = zeros(n,n);
    
    % array holding variance (for convergence detection)
    v = zeros(n,n);
    m = zeros(n,n);
    % variance from last round
    old_v = v;

    
    %iterate
   for k=1:max_iter

      %variable to check nodes
      % for each check
      for i=1:n
         % for each var
           for j=1:n
               %if there is an edge from var node i to check node j
               if ((H(i,j) ~= 0) && (i~=j))
                   
                  if (k == 1)
                      %v(j,i) = -H(j,i).^2/sigma;
                      %M{j,i} = kde(-H(j,i)*y(i)/v(j,i), sqrt(H(j,i).^2/sigma), 1);% 
                      prod = self_pot{i};
                      verify_pdf(prod);
                      %prod = shiftpdf(prod,y(i),xticks);
                      
                      
                      M{i,j} = pdf_integral(prod, model_order, -H(i,j), xticks, delta, epsilon);
                      verify_pdf(M{i,j});

%                       maxval = find(M{i,j} == max(M{i,j}));
%                       if (length(maxval) > 1)
%                           m(i,j) = mean(xticks(maxval));
%                       else
%                           m(i,j) = xticks(maxval);
%                       end
                      %M2N{in,neighbors(j)} = pdf_integral(pdf_res, model_order, A(in,neighbors(j)),xticks,delta,epsilon);
                      %mji = -H(j,i)*y(i)/sigma;
                      m(i,j) = y(i)*-H(i,j);
                      M_shift(i,j) = y(i)*-H(i,j);
                  else      
                      % node degree
                      dd = sum(H(:,i)~= 0);
                      assert(H(i,i) == 0);
                      toMul = cell(1, dd);
                      toMul{1} = self_pot{i};% 
                      toAdd = 0;
                      
                      cc = 2;
                      for l=1:n
                          if ((H(l,i) ~= 0) && (l~=j))
                              current = old_M{l,i};
                              verify_pdf(current);
                              toMul{cc} = current;
                              toAdd = toAdd + old_M_shift(l,i);
                              cc = cc + 1;
                          end
                      end % for l

                    % computes the approx. product of incoming mixtures
                     prod = mulpdfc(toMul);
                     verify_pdf(prod);
                     %prod = shiftpdf(prod,y(i),xticks);
                     
                     %prod = prod/sum(prod);


                     
                     prod2 = mulpdf(prod, old_M{j,i});
                     [mrecon(i), srecon(i),xrecon(i)]=meanvarmaxpdf(prod2, xticks); % computes statistics   
                     mrecon(i) = y(i) + toAdd + old_M_shift(j,i);
                     %
                      %M{i,j} = integral_ekde(prod, H(i,j));
                      M{i,j} = pdf_integral(prod, model_order, -H(i,j), xticks, delta, epsilon);
                      verify_pdf(M{i,j});
                      m(i,j) = -H(i,j)*(y(i) + toAdd);
                      M_shift(i,j) = -H(i,j)*(y(i) + toAdd);
                      disp([num2str(i) '.' num2str(j) ' hprod ' num2str((y(i) + toAdd))]);
                  end
                  
                  
               end % if
           end % for j
          
      end % for i

      m
      old_M = M;
      old_M_shift = M_shift;
      
      if (k >= 2)
          dispvec_anderrors_gabp(mrecon, H, y, dispind, x);
      end
      %disp(['NBP iteration ', num2str(k)]);      
   end
           
  %variable to check nodes
      % for each check
      xx=cell(1,n);
      mm=zeros(1,3);
      vv=mm;
      mx=mm;
      for i=1:n
         % for each var
               %if there is an edge from var node i to check node j
                
                 % node degree
                  dd = sum(H(:,i)~= 0);
                  toMul = cell(1, dd+1);
                  toAdd = 0;
                  toAdd2 = 0;
                  cc = 1;
                  for l=1:n
                      if ((H(l,i) ~= 0) && (l~=i))
                          current = M{l,i};
                          verify_pdf(current);
                          %current = reverse_gabp(current,model_order);
                          toMul{cc} = current;
                          toAdd = toAdd + xticks(current == max(current));
                          toAdd2 = toAdd2 + sum(xticks.*current);
                          cc = cc + 1;
                      end
                  end % for l
                  toMul{cc} = self_pot{i};% 
                  toAdd = toAdd + y(i);
                  toAdd2 = toAdd2 + y(i);
                  
                  % computes the approx. product of incoming mixtures
                  prod = mulpdfc(toMul);
                  verify_pdf(prod);
                  prod = shiftpdf(prod,y(i),xticks);
                  [meanx(i),varx(i),maxx(i)]=meanvarmaxpdf(prod,xticks);
                  meanx(i) = toAdd;
                  maxx(i) = toAdd2;
      end % if
      
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
%        new_mean = sum(xticks.*npdf);
%        if (length(new_mean) == 1 && (abs(new_mean - old_mean * edge) >= 2*delta))
%            %assert(abs(new_mean - old_mean * edge) < delta);%DB TODO 
%        end
    end
    verify_pdf(npdf);
end
  
% sanity checks
function []=verify_pdf(pdf)
   assert(sum(pdf < 0)==0);
   assert(sum(isnan(pdf))==0);
   assert(sum(isinf(pdf))==0);
   assert(length(unique(pdf)) > 1);
   assert(sum(~isreal(pdf)) == 0);
end


% point by point multiplication
function c=mulpdf(a,b)
c=a.*b;
%c=c/sum(abs(c));
end

function [ret]=mulpdfc(a)
    len=length(a);
    ret = a{1};
    verify_pdf(ret);
    for pp=2:len
       verify_pdf(a{pp});
       ret = mulpdf(ret,a{pp});
       verify_pdf(ret);
    end
    ret=ret/sum(ret);
end

% function c=divpdf(a,b)
% b=b+1e-100;
% c=a./b;
% %c=c/sum(abs(c));
% end
% 
% function [ret]=divpdfc(a)
%     len=length(a);
%     ret = a{1};
%     verify_pdf(ret);
%     for pp=2:len
%        verify_pdf(a{pp});
%        ret = divpdf(ret,a{pp});
%        verify_pdf(ret);
%     end
%     ret=ret/sum(ret);
% end


function [m,sig,mp]=meanvarmaxpdf(pdf, xx)
    pdf=pdf/sum(pdf);
    [mm,ind]=max(pdf);
    mp=xx(ind);
    m=sum(xx.*pdf);
    sig=sqrt(sum(xx.*xx.*pdf)-m*m);
end


function [] =dispvec_anderrors_gabp (v, A, y, dispind, x)
    l=length(dispind);
    nn=length(y);
    s='[';
    for pp=1:l
      s=sprintf('%s %7.2f',s,v(dispind(pp)));
    end
    s=sprintf('%s]',s);
    %mv=A(1:nn,nn+1:2*nn)*v(1:nn)';
    mv=A*v';
    er=norm(y-mv(1:nn));
    %tnorm = norm(x'-v(1:length(y)));
    tnorm = norm(x'-v);
    s=sprintf('%s (%7.4f) (%7.4f)',s, er, tnorm);
    disp(s);
end

% DB: Since the array is fftshifted, the first position is zero (odd sample size assumed), no need to
% reverse it. The other values are reversed, which means we computed
% integral with edge weight of -1.
function op=reverse_gabp(ip, m)
    op=[ip(m:-1:1)];
end


%move the new mean to "Offset - old mean"
function pdf=shiftpdf(pdf, offset, xx)
    verify_pdf(pdf);
    pdf = interp1(xx,pdf,xx-offset);
    pdf(isnan(pdf)) = 0;
    pdf = pdf./sum(pdf);
end

end

