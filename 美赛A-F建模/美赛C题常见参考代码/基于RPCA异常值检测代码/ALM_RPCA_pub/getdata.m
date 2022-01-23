function opts = getdata(dataformat)

if strcmp(dataformat,'surveillance-video-Hall')
    load Hall_airport_1000_1496_497_144_176_gray.mat;
    D = images(:,1:200); [n1,n2] = size(images); % imn1 = 144; imn2 = 176;
    opts.D = D; opts.mu = norm(D)/1.25;
    Xs = D; Ys = D;
end

opts.Xs = Xs;  opts.Ys = Ys;
opts.n1 = n1; opts.n2 = n2;
opts.sigma = 1e-6; opts.maxitr = 500; opts.rho = 1/sqrt(n1); 
opts.eta_mu = 2/3; opts.eta_sigma = 2/3;
opts.muf = 1e-6;
opts.sigmaf = 1e-6;
opts.epsilon = 1e-7;
opts.sv = 100;
