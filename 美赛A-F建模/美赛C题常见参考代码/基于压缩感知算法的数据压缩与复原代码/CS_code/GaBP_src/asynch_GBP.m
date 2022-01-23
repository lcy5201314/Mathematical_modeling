%This is an implementation of the Gaussian BP algorithm (async version)
% Written by Danny Bickson
%See: http://books.nips.cc/papers/files/nips18/NIPS2005_0210.pdf
%Equations 7,8,9
%Input: A - information matrix mxm, (assumed to be symmetric) and
% diagonally dominant.
% b - shift vector 1xm
% maxround - maximum round number
% epsilon - threshold for convergence 
% quiet (optional) - if true does not print to screen
%Output: The solution for the inference problem
%        vector h of size 1xm s.t. h = max(1/2h'Ah +h'b)
%        J - vector of the values Pii (the diagonal of the matrix A^-1)
%        r - round number of convergence (if converged)

function [h,J,r,C] = asynch_GBP(A,b,maxround,epsilon,quiet)
m=length(A);
C=zeros(maxround,length(b));
%messages
Mh=zeros(m,m);
MJ=zeros(m,m);

%return values
h=zeros(1,m);
J=zeros(1,m);

conv = false;

if ~exist('quiet','var')
   quiet = false; 
end

% algorithm rounds
for r=1:maxround
    %disp(['starting async GBP round ', num2str(r)]); 
    old_MJ = MJ;
    old_Mh = Mh;
   
	% for each node
	for i=1:m
		% sum up all mean and percision values got from neighbors
		h(i) = b(i) + sum(Mh(:,i));  %(7)
	    
        %variance can not be zero (must be a diagonally dominant matrix)!
        %assert(A(i,i) ~= 0);
        J(i) = A(i,i) + sum(MJ(:,i));
	
		% send message to all neighbors
		for j=1:m
			if (i ~= j && A(i,j) ~= 0)
				h_j = h(i) - Mh(j,i);
				J_j = J(i) - MJ(j,i);

                %assert(A(i,j) == A(j,i));
			    %assert(J_j ~= 0);
                Mh(i,j) = (-A(j,i) / J_j)* h_j; %(8)
				MJ(i,j) = (-A(j,i) / J_j) * A(i,j);
            end
        end
        
        %h
    end

    
    for i=1:m
        h(i) = b(i) + sum(Mh(:,i));  %(9)
        J(i) = A(i,i) + sum(MJ(:,i));
    end
    J=1./J;
    h=h.*J;
    
    C(r,:)=h;
        
    
    max_diff  = epsilon;
   
    if (r > 2 && ((norm(C(r,:) - C(r-1,:))/norm(C(r,:))) < epsilon))
        if ~quiet
        disp(['Async GBP (MJ) Converged afeter ', num2str(r), ' rounds ']); 
        end
        
        for i=1:m
            h(i) = b(i) + sum(Mh(:,i));  %(9)
            J(i) = A(i,i) + sum(MJ(:,i));
        end
       conv = true;
        break;
    end

end



if (conv == false)
    disp(['GBP (MJ) did not converge in ', num2str(r), ' rounds ']); 
    for i=1:m
        h(i) = b(i) + sum(Mh(:,i));  %(9)
        J(i) = A(i,i) + sum(MJ(:,i));
    end

end




J = 1./J;
h=h.*J;
end




