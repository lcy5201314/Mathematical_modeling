[Y,fs,bits]=wavread('wrd.wav');
% plot(Y);
% hold on;

M1=30:40:280;
M2=M1-20;

for ii=1:7
    ii
    [r(ii),snr(ii)]=decisionomp(Y,M1(ii),M2(ii));
end
hold off
figure(2)
plot(r,snr,'-o')
title('decision omp')