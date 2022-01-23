%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a sample code for testing the algorithm in our paper 
% ''Fast Alternating Linearization Methods for
%       Minimizing the Sum of Two Convex Functions'', Donald Goldfarb,
%       Shiqian Ma and Katya Scheinberg, Tech. Report, Columbia University,
%       2009 - 2010. Preprint available at: 
%       http://arxiv.org/pdf/0912.4571v2.pdf
%
% Author: Shiqian Ma
% Date  : Apr. 20, 2010 
% IEOR, Columbia University, Copyright (2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get data 
tic;
dataformat = 'surveillance-video-Hall';

opts = getdata(dataformat); 
time_getdata = toc;
fprintf('%f seconds to get data ! \n', time_getdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Call ALM to solve the problem
tic; out_ALM = ALM_SADAL_smoothed(opts.D,opts); time_ALM = toc;
fprintf('*******************************************************************\n');
%%%%%%%%% print stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('*******************************************************************\n');
fprintf('ALM  : iter: %d, StopCrit: %3.2e, time: %f\n', ...
    out_ALM.iter, out_ALM.StopCrit, time_ALM);


% plot the images 
D = opts.D; X = out_ALM.X; Y = out_ALM.Y; 
imn1 = 144; imn2 = 176; ind1 = 30; ind2 = 96; ind3 = 159;
subplot(3,3,1); imshow(reshape(D(:,ind1),imn1,imn2),[]); axis off; title('Video')
subplot(3,3,2); imshow(reshape(X(:,ind1),imn1,imn2),[]); axis off; title('Background')
subplot(3,3,3); imshow((reshape(Y(:,ind1),imn1,imn2)),[]); axis off; title('Foreground')
subplot(3,3,4); imshow(reshape(D(:,ind2),imn1,imn2),[]); axis off; 
subplot(3,3,5); imshow(reshape(X(:,ind2),imn1,imn2),[]); axis off; 
subplot(3,3,6); imshow((reshape(Y(:,ind2),imn1,imn2)),[]); axis off; 
subplot(3,3,7); imshow(reshape(D(:,ind3),imn1,imn2),[]); axis off; 
subplot(3,3,8); imshow(reshape(X(:,ind3),imn1,imn2),[]); axis off; 
subplot(3,3,9); imshow((reshape(Y(:,ind3),imn1,imn2)),[]); axis off; 
