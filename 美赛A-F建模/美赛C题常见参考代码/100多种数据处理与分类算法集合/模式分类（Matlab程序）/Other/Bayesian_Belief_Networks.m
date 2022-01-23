function [decision, P] = Bayesian_Belief_Networks(net, data)

% Find the most likely decision given a Bayesian belief network and data for the decision
%
% Inputs:
%	net     - A bayesian belief network (See below for the structure)
%   data    - The known data from which to estimate
%
% Output:
%	decision- The most probable decision based on this data
%   P       - Probabilites from which the decision was obtained
%
% The network should be a vector of structures such that each cell in the net vector is built as follows:
% net(i).child  = <Child node number/s (Must be different from i), or [] if this is a termination node
% net(i).P      = <Node probability matrix given the parents, i.e. P(i|j,k)>
% net(i).comment= <Optional comment>
%
% Data should be a 1xN vector, where N is the number of nodes in the tree
% This vector should contain a number where the choice is known, 0 where it is unknown and inf 
% where the decision is to be taken
%
% NOTE: This implementation takes a maximum of two parents and children per node.
%
% See also the net provided in the file bayes_belief_net, which enumerates example 4 in DHS chapter 2

Nnodes  = length(net);
if (Nnodes ~= size(data,2)),
    error('Mismatch between data length and network size')
end

%First, find who the parent nodes are and invert the net, that is, for each node, locate it's parents
%A parent node is a node that is not pointed to by other nodes
parents = [];
dims    = zeros(1,Nnodes);  %This will hold the number of values for each node
for i = 1:Nnodes,
    net(i).parents  = [];
    dims(i)         = size(net(i).P,1);
    for j = 1:Nnodes,
        if ~isempty(net(j).child),
            child_j = find(net(j).child == i);
            if ~isempty(child_j),
                net(i).parents = [net(i).parents, j];
            end
        end
    end
    
    if isempty(net(i).parents),
        parents = [parents, i];
    end
end
disp(['The parent nodes are: ' num2str(parents)])

%Introduce the evidence into the network
for i = 1:Nnodes,
    if ((data(i) < inf) & (data(i) > 0)),
        indices = find(~ismember(1:size(net(i).P,1),data(i)));
        net(i).P(indices,:,:) = 0;
        net(i).P(data(i),:,:) = net(i).P(data(i),:,:) ./ sum(sum(squeeze(net(i).P(data(i),:,:))));
    end

    %Introduce two variables which will be needed later
    net(i).Pp       = inf;
    net(i).Pc       = inf;
end

%Now, find the probability of each choice for each node in the network
up  = ones(1, Nnodes);
down= ones(1, Nnodes);
target = find(~isfinite(data));

while ((sum(up) + sum(down) > 0) & ((up(target)==1) | (down(target)==1))),
   for i = 1:Nnodes,
      %Can we stop?
      if (up(target)==0) & (down(target)==0),
         break
      end
      
        if up(i),
            %Are the parents computed?
            if ~isempty(net(i).parents),
               if sum(up(net(i).parents)~=0)==0,
                  nodeP = net(i).P;
                    %The parents are computed, so compute this node's probability
                    switch length(net(i).parents),
                    case 1,
                        net(i).Pp= sum((ones(size(nodeP,1),1)*net(net(i).parents).Pp).* nodeP,2)';
                    case 2,
                        temp		= net(net(i).parents(1)).Pp' * net(net(i).parents(2)).Pp;
                        for j = 1:size(net(i).P,1),
                           net(i).Pp(j) = sum(sum(squeeze(nodeP(j,:,:)).*temp,1),2)';
                        end
                    otherwise,
                    end
                    up(i) = 0;
                end
            else
                net(i).Pp   = net(i).P';
                up(i)       = 0;
            end
        end
        
        if down(i),
            %Are the childrens computed?
            if ~isempty(net(i).child),
                if sum(down(net(i).child)~=0)==0,
                    %The children are computed, so compute this node's probability
                    nodeP = net(i).P;
                    switch length(net(i).child),
                    case 1,
                        net(i).Pc= sum(net(net(i).child).Pc);
                    case 2,
                        net(i).Pc= sum(net(net(i).child(1)).Pc) .* sum(net(net(i).child(2)).Pc);
                    otherwise,
                    end
                    down(i) = 0;
                end
            else
                net(i).Pc   = net(i).P;
                down(i)     = 0;
            end
        end
    end
end


%Finally, compute the node probability
target  = find(~isfinite(data));
P       = net(target).Pp .* net(target).Pc;
P       = P ./ sum(P);
decision= find(P == max(P));