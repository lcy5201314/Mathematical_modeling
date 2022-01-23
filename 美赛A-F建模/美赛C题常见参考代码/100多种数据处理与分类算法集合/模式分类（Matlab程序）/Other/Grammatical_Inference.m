function [A,I,S,P] = Grammatical_Inference(x, labels)

% Bottom-Up Parsing
%
% Inputs:
%	x					- Text vectors to parse
%   labels              - Labels for the text vectors {-1, 1}
%
% Output:
%	A					- Alphabet vector
%	I					- Variables vector
%	S					- Root symbol
%	P					- Production rules
%
% NOTE: THis is a very simple implementation. You may 
% want to change the rule-adding methods to obtain better results according to your data


%Define the root symbol
S       = 'S';

%First, identify the alphabet vector (Very easy)
X       = char(x(find(labels)));
A       = unique(X(:));
A       = A(find(A~=32));

%The variable vector will be defined simply as upper case of the alphabet vector
I       = upper(A);

%And now to the difficult part: Finding the production rules
P       = cellstr(['S->' I(1)]);
%P(2)    = {'S->S'};
n_plus  = sum(labels>0);
in_plus = find(labels ==  1);
in_minus= find(labels == -1);

for i = 1:n_plus,
    %Read x_i_plus from Dplus
    data        = char(x(in_plus(i)));
    plus_parsed = 0;
    
    while (~plus_parsed),
          
        tempP = {};
        
        %If x_i_plus cannot be parsed by G
        [V, success] = Bottom_Up_Parsing(A,I,S,P,data);
        if success,
            plus_parsed = 1;
            break
        end
        
        if isempty(V(1).string),
            %Propose additional productions to P:
            %1. If V is empty, add simple grammar rules so that every letter inthe alphabet is converted to a variable
            found   = 0;
            for j = 1:length(A),
                for k = 1:length(P),
                    if ~isempty(findstr(A(j), char(P(k))))
                        found = 1;
                        break
                    end
                end
                if ~found,
                    tempP(length(tempP)+1) = cellstr([I(j) '->' A(j)]);
                end
            end
        else
            %Propose additional productions to P:    
            %2. If V is not empty, add more complicated rules
            %We do it by adding simple rules that connect the last layers of A to the output
            Np      = length(P);
            
            %Find an empty space to merge into
            found = 0;
            for j = 2:length(data),
                slot = 0;
                for i = 1:length(data)-j+1,
                    if isempty(V(i,j).string),
                        slot = j;
                        line = i;
                        found= 1;
                        break
                    end
                end
                if found,
                    break
                end
            end
            
            %Find who to merge that doesn't have a rule already
            for i = 1:slot-1,
                B   = V(line, i).string;
                C   = V(i+line,slot-i).string;
                rule= strcat(B, C);
                if (rand(1) > .5),
                    target = B;
                else
                    target = C;
                end
                rule = cellstr([char(target) '->' char(rule)]);
                
                %Find that this rule doesn't exist yet
                found = 0;
                for k = 1:length(P),
                    if strmatch(char(P(k)), char(rule)),
                        found = 1;
                    end
                end
                
                if ~found,
                    tempP(length(tempP)+1) = rule;
                end
            end
            
        end

        %Accept updates if G parses x_i_plus but no string in D_minus
        parse = 0;
        for j = 1:length(in_minus),
            [V_minus, success] = Bottom_Up_Parsing(A,I,S,[P, tempP],x(in_minus(j)));
            if success
                parse = 1;
                break
            end
        end
        
        %If not parsed, add the new rules
        if ~parse,
            P = [P, tempP];
        end
        
        %Check if the new rule has solved the problem
        [V, success] = Bottom_Up_Parsing(A,I,S,P,data);
        if success,
           plus_parsed = 1;
        end


    end
end