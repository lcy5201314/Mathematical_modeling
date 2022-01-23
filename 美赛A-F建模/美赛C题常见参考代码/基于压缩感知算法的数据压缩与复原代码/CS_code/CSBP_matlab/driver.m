%-------------------
% initializations
%-------------------
epsilon=1e-70;
aux_rows=l; 
disp (sprintf('n = %d    k = %d   l = %d', n, k, l));

%-------------------
% create auxiliary data structures
%-------------------
[aux, aux_rows_actual]=get_aux(phi, phisign, n, l, aux_rows);
[self_indexN, self_indexM]=GetSelfIndices(phi, aux);   %Used in BP
[tmptmp, dispind]=sort(-abs(x));
dispind=[dispind(1:5), 1:5];
  
%-------------------
% Encode (compute measurements)
%-------------------
disp ('GENERATING THE MEASUREMENTS...');
measvec=encoder(phi, phisign, x);
disp (sprintf('Squared l2 norm of x : %g', (norm(x))^2)); 
disp (sprintf('Squared l2 norm of y_0 : %g', (norm(measvec))^2) ); 
disp (sprintf('msmt snr set: %g', SNR)); 
disp (sprintf('msmt snr actual: %g', ((norm(measvec))^2)/(length(measvec)) )); 
disp(sprintf('  ... Number of measurements=%d', length(measvec)));
measvec=measvec+ynoise_sigma_actual*randn(length(measvec),1);	% add noise

%-------------------
% Decode
%-------------------
disp ('STARTING THE DECODER...');
boundx=sigma_signal*10;

ynoise_sigma=ynoise_sigma_actual;
[xrecon, mrecon, srecon, pdf, pdf_xx, pdf_prior]=decoder(x, measvec, n, k, l, max(aux_rows_actual), sigma_signal, sigma_noise,ynoise_sigma,iter, phi, phisign, aux, model_order, self_indexN, self_indexM, dispind, boundx, epsilon, gamma_pdbp, gamma_mdbpf, gamma_mdbpb);
