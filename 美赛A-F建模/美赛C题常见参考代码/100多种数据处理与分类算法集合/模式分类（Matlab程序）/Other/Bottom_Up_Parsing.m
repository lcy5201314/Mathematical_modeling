function [V, success] = Bottom_Up_Parsing(A,I,S,P,x)

% Bottom-Up Parsing
%
% Inputs:
%	A					- Alphabet vector
%	I					- Variables vector
%	S					- Root symbol
%	P					- Production rules (NOTE: Do not merge operators using the OR operator)
%	x					- Text vector to parse
%
% Output:
%	V					- Parsing table. 

%NOTE: Because of the matlab indexing, which starts from 1, the algorithm is slightly
%different from the pseudocode in DHS

%Translate the rules matrix into rules and targets
for i = 1:length(P),
   t			= char(P(i));
   l			= findstr(t, '->');
   if isempty(l),
      error('Incorrect rule. A rule must contain the string ''->''')
   else
      Ptargets(i) = cellstr(t(1:l-1));
      P(i)		 = cellstr(t(l+2:end));
   end
end

n	= length(x);
for i = n:-1:1,
   for j = n:-1:1,
      V(i,j).string = '';
   end
end

for i = 1:n,
   %V(i,1) <- {A|A->x(i)}
   rule				= strmatch(x(i), P, 'exact');
   V(i,1).string	= Ptargets(rule);
end

%Do another pass for cases such that P={'S->A','A->a'}
if n == 1,       
	for i = 1:length(P),
        if strmatch(Ptargets(i),'S'),
            if strmatch(P(i), V(j,1).string),
                V(j,1).string = 'S';
            end
        end
    end
end

for j = 2:n,
   for i = 1:n-j+1,
      V(i,j).string = '';
      for k = 1:j-1,
         %V(i,j)<-V(i,j)U{A|A->BC in P, B in V(i,k), C in V(i+k,j-k)
         B	= V(i,k).string;
         C	= V(i+k, j-k).string;
         for b = 1:length(B),
            for c = 1:length(C),
                rule 	= strcat(B(b), C(c));
                target= strmatch(rule, P, 'exact');
                if ~isempty(target),
                   V(i,j).string = [V(i,j).string, Ptargets(target)];
                end
            end
         end
      end
   end
end

success = 0;
if ~isempty(strmatch(S,V(1,n).string,'exact')),
    if nargout == 2,
        success = 1;
    else
        disp(['Parse of ' x ' successful in G'])
    end
end
