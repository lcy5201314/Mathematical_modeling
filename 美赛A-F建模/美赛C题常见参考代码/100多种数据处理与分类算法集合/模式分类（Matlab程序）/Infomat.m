function gamma = infomat(patterns, targets)

%Calculate the mutual information matrix (Use by the Koller algorithm)
%
%Inputs:
%	patterns - Input patterns
%   targets  - Input targets (0 or 1)
%
%Outputs:
%	gamma    - The information matrix

[Nf, L]	    = size(patterns);
Kdist		= 2;				%How many patterns to group together
Nhist		= floor(L^(1/3));	%How many bins for the histogram
Uc          = unique(targets);

%Normalise patterns to [-1, 1]
patterns    = (2*patterns - (min(patterns')+max(patterns'))'*ones(1,L)) ./ ((max(patterns')-min(patterns'))'*ones(1,L));

gamma		= zeros(Nf);
P           = zeros(1, length(Uc));
for i = 1:length(Uc),
    P(i) = length(find(targets == Uc(i)));
end

disp('Started calculation of the cross-entropy matrix')

for i = 1:Nf,
   tic
   hist_i   = high_histogram(patterns(i,:),Nhist);
   for j = 1:length(Uc)
       in             = find(targets == Uc(j));
       g_j_d_i(j,:,:) = (high_histogram(patterns(i,in),Nhist).*hist_i)*ones(1,Nhist);
   end

   for j = i:Nf,
       hist_ij= high_histogram(patterns([i,j],:),Nhist);
       hist_j = high_histogram(patterns(j,:),Nhist);
       for k = 1:length(Uc),
           in              = find(targets == Uc(k));
           g_k_n(k,:,:)    = high_histogram(patterns([i,j],in),Nhist).*hist_ij;
           g_k_d_j(k,:,:)  = (high_histogram(patterns(j,in),Nhist).*hist_j)*ones(1,Nhist);
       end 
      
      %The addition of eps and the multipication by sign is for numeric reasons
      for k = 1:length(Uc),
          g_k_j(k,:,:) = squeeze(g_k_n(k,:,:))/P(k).*log(eps + squeeze(g_k_n(k,:,:))./...
                        (squeeze(g_k_d_j(k,:,:)) + eps)).*sign(squeeze(g_k_d_j(k,:,:)));
          g_k_i(k,:,:) = squeeze(g_k_n(k,:,:))/P(k).*log(eps + squeeze(g_k_n(k,:,:))./...
                        (squeeze(g_j_d_i(k,:,:)) + eps)).*sign(squeeze(g_j_d_i(k,:,:)));
      end
      gamma(i,j) = squeeze(sum(sum(sum(g_k_j))));
      
      gamma(j,i) = squeeze(sum(sum(sum(g_k_i))));
   end
   
   t = toc;  
   disp(['Iteration ' num2str(i) ': Time taken: ' num2str(t) '[sec]'])
end
