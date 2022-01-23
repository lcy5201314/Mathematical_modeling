function delta = CARTfunctions(split_point, patterns, targets, dim, split_type)

%Calculate the difference in impurity for the CART algorithm
Uc	= unique(targets);
for i = 1:length(Uc),
   in		= find(targets == Uc(i));
   Pr(i)	= length(find(patterns(dim, in) >   split_point))/length(in);
   Pl(i)	= length(find(patterns(dim, in) <=  split_point))/length(in);
end

switch split_type,
case 'Entropy'
	Er		= sum(-Pr.*log(Pr+eps)/log(2));
   El		= sum(-Pl.*log(Pl+eps)/log(2));
case {'Variance', 'Gini'}
	Er		= 1 - sum(Pr.^2);
   El		= 1 - sum(Pl.^2);
case 'Missclassification'
   Er		= 1 - max(Pr);
   El		= 1 - max(Pl);
otherwise
   error('Possible splitting rules are: Entropy, Variance, Gini, or Missclassification')
end

P		= length(find(patterns(dim, :) <= split_point)) / length(targets);
delta = -P*El - (1-P)*Er;
