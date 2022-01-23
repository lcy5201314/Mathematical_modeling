function [test_indices, train_indices] = make_a_draw(how_many, out_of_how_many)

%Choose a number of indices out of the maximum number allowed.

indices		  = randperm(out_of_how_many);
train_indices = sort(indices(1:how_many));
test_indices  = sort(indices(how_many+1:out_of_how_many));
