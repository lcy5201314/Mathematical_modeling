function mrecon=driver(n,k,l,phi,phisign,x,SNR,sigma_1,sigma_0,sigma_Z,...
   iter,p,gamma_mdbpf,gamma_mdbpb,gamma_pdbp)
%-------------------
% driver_function.m
% This is a function version of the script file driver.m
% Inputs:
% n - length of signal
% K - number of large coefficients
% l - number of ones per row in LDPC-CS matrix
% phi - encoding matrix structure
% phisign - signs of nonzero matrix entries
% x - actual signal
% SNR - input SNR
% sigma_1 - large signal coefficients
% sigma_0 - small coefficients
% sigma_Z - measurement noise
% iter - number of BP iterations
% p - number of samples
% gamma_mdbpf,gamma_mdbpb,gamma_pdbp - parameters used for damping
% Output: mrecon - estimated value of xhat.
% Dror, 12.22.2008
%-------------------

%-------------------
% create auxiliary data structures
%-------------------
aux_rows=l; 
[aux, aux_rows_actual]=get_aux(phi, phisign, n, l, aux_rows);
[self_indexN,self_indexM]=GetSelfIndices(phi, aux);   %Used in BP
[tmptmp, dispind]=sort(-abs(x));
dispind=[dispind(1:5), 1:5];
  
%-------------------
% Encode (compute measurements)
%-------------------
disp ('GENERATING THE MEASUREMENTS...');
measvec=encoder(phi, phisign, x);
disp(sprintf('Number of measurements=%d', length(measvec)));
measvec=measvec+sigma_Z*randn(length(measvec),1);	% add noise

%-------------------
% Decode
%-------------------
disp ('STARTING THE DECODER...');
epsilon=1e-70;
boundx=sigma_1*10;
[xrecon, mrecon, srecon, pdf, pdf_xx, pdf_prior]=decoder(x,measvec,n,k,l,max(aux_rows_actual),sigma_1,sigma_0,sigma_Z, iter, phi, phisign, aux,p,self_indexN, self_indexM, dispind, boundx, epsilon, gamma_pdbp, gamma_mdbpf, gamma_mdbpb);
