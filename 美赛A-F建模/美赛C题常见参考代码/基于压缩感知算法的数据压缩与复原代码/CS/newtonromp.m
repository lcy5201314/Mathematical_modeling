 function [r,vOut, numIts] = newtonromp(n, Phi, y,lamada)
% [vOut] = romp(n, Phi, y)
%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% d = ambient dimension of the signal x
% N = number of measurements
% n = sparsity level of n
% Phi = N by d measurement matrix
%y = measurement vector (Phi * x)
% vOut = reconstructed signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION DESCRIPTION %%%%%%%%%%%%%%%%%
% romp takes parameters as described
% above. Given the sparsity level n and
% the N by d measurement matrix Phi, and
% the measurement vector y = Phi * x, romp
% reconstructs the original signal x.
% This reconstruction is the output.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear r I J J0 u b ix numIts Jvals
warning off all
N = size(Phi, 1);
d = size(Phi, 2);
% Set residual
r = y;
%Set index set to "empty"
I = [];xs=zeros(d,1);
%Counter (to be used optionally)
numIts = 0;

%Run ROMP
while length(I)-1 < 2*n && norm(r) > 10^(-3)
    numIts = numIts + 1;
    %Find J, the biggest n coordinates of u
    u = Phi' * r;
    absu = abs(u);
    [b, ix] = sort(absu, 'descend');
    J = ix(1:n);
    Jvals = b(1:n);
    %Find J0, the set of comparable coordinates with maximal energy
    w=1;
    best = -1;
    J0 = zeros(1);
    average_energy=0;
    while w <= n
        first = Jvals(w);
        firstw = w;
        energy = 0;
        while ( w <= n ) && ( Jvals(w) >=lamada* first )
            energy = energy + Jvals(w)^2;
            w = w+1;
        end
%         if energy > best
%             best = energy;
%             J0 = J(firstw:w-1);
%         end
        av_e=energy/length(first:w-1);
        if energy>best
        if  2*av_e>average_energy
            
            J0 = J(firstw:w-1);
            best = energy;
            average_energy=av_e;
            
        end
        end
    end
    %Add J0 to the index set I
    I(length(I)+1: length(I)+length(J0)) = J0;
   
    %Update the residual
    PhiSubI = Phi(:, I(1));
    for c=2:length(I)
        if ~isMember(I(1:c-1),I(c))
            PhiSubI(:,c) = Phi(:, I(c));
        end
    end
%     Aug_t=PhiSubI;
%     gn=Aug_t'*r;   %i*1,梯度
%     cn=Aug_t*gn;    %M*1
%     delta=(r'*cn)/norm(cn).^2;%步长
%     r=r-Aug_t*delta*gn; 
   Aug_t=PhiSubI;
    H=-Aug_t'*Aug_t;                              %海森矩阵i*i
    gn=Aug_t'*r;                                   %i*1,梯度
   delta=-H^(-1)*gn;                                      %方向
 a=1;%步长
   xs(I)=xs(I)+a*delta;
 r=r-a*Aug_t*delta;   

%  r=r-Aug_t*xs;
 vOut = xs;
% error=norm(hat_x-x)^2/norm(x)^2  ;
% snr=10*log10(1/error)
%     x= lscov(PhiSubI, y);
%     r = y - PhiSubI * x;

end 
I_length=length(I)

    

