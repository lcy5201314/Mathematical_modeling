function [r,vOut, numIts] = romp(n, Phi, x)
% [vOut] = romp(n, Phi, x)
%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% n = sparsity level of n




vOut = zeros(size(Phi, 2), 1);
% Set residual
r = x;
%Set index set to "empty"
I = [];
%Counter (to be used optionally)
numIts = 0;
%Run ROMP
while length(I)-1 <2*n && norm(r) > 10^(-6)
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
        while ( w <= n ) && ( Jvals(w) >= 0.5* first )
            energy = energy + Jvals(w)^2;
            w = w+1;
        end
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
    I=unique(I);
    %Update the residual
    PhiSubI = Phi(:, I);
    
    y = lscov(PhiSubI, x);
    r = x - PhiSubI * y;
end % end Run IRA
vSmall = PhiSubI \ x;

vOut(I)=vSmall;
Iromp=length(I);
