function test_targets = Genetic_Algorithm(train_patterns, train_targets, test_patterns, params)

% Classify using a basic genetic algorithm
% Inputs:
% 	training_patterns   - Train patterns
%	training_targets	- Train targets
%   test_patterns       - Test  patterns
%	Params              - [Type, TargetErr, Nchrome, Pco, Pmut], where:
%                           Type:       Type of weak learner 
%                           Params:     Parameters of weak learner 
%                           TargetErr:  Target error on the train set for the GA
%                           Nchrome:    Number of chromosomes to use
%                           Pco:        Probability of recombination
%                           Pmut:       Probability of mutation
%
% Outputs
%	test_targets        - Predicted targets

[type, Wparams, TargetErr, Nchrome, Pco, Pmut] = process_params(params);
[D,L]     = size(train_patterns);
iter      = 0;

%Build the chromosomes
%The mapping in this realization is wheather or not to use a given example for building the classifier
chromosomes = rand(Nchrome, L)>0.5;
ranking     = ones(1,Nchrome);

while 1,    
    %Determine the fit of each chromosome
    for i = 1:Nchrome,
        if (ranking(i) == 1),
            %Build a classifier and test it
            index       = find(chromosomes(i,:) == 1);
            targets     = feval(type, train_patterns(:,index), train_targets(index), train_patterns, Wparams);
            ranking(i)  = mean(xor(targets, train_targets));
        end
    end
    
    if ((min(ranking) < TargetErr)),
        break
    end

    iter = iter + 1;
    if (iter/10 == floor(iter/10)),
        disp(['Iteration number ' num2stR(iter) ': Best classifier so far has an error of ' num2str(min(ranking)*100) '%'])
    end

    %Rank the chromosomes
    [m, rating] = sort(ranking);
    
    for i = 1:floor(Nchrome/2),
        %Select the two chromosomes with the highest score
        c1  = rating(i*2-1);
        c2  = rating(i*2);
        
        %If rand[0,1]<Pco then 
        if (rand(1) < Pco),
            %Crossover each pair at a random bit
            crossover = randperm(L-2);  %This is to avoid edges
            temp1     = chromosomes(c1,:);
            temp2     = chromosomes(c2,:);
            chromosomes(c1,1:crossover(1)+1)    = temp2(1:crossover(1)+1);
            chromosomes(c2,crossover(1)+2:end)  = temp1(crossover(1)+2:end);
            ranking(c1) = 1;
            ranking(c2) = 1;
        else
            %Change each bit with probability Pmut
            chromosomes(c1,:) = xor(chromosomes(c1,:),(rand(1,L)>(1-Pmut)));
            chromosomes(c2,:) = xor(chromosomes(c2,:),(rand(1,L)>(1-Pmut)));
            ranking(c1) = 1;
            ranking(c2) = 1;
        end
    end
    
end

best            = find(ranking == min(ranking));
best            = best(1);
index           = find(chromosomes(best,:) == 1);
test_targets    = feval(type, train_patterns(:,index), train_targets(index), test_patterns, Wparams);