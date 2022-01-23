function indices = combinations (in_from, N)

% Find all N combinations of input indices 
%
% Inputs:
%   in_from - The input indices
%   N       - How many indices are requiered
%
% Outputs:
%   indices - The indice combinations

Nf = length(in_from);

comb_mat  = zeros(2^Nf, Nf);
basic_mat = [0, 1];

for i = 1:Nf,
    rep_mat     = ones(2^(i-1), 1) * basic_mat;
    pattern     = rep_mat(:);
    rep_pattern = pattern * ones(1, 2^(Nf-i));
    column      = rep_pattern(:);
    comb_mat(:,Nf - i + 1) = column;
end

%Find the number of the right rows
in = find(sum(comb_mat') == N);

indices = zeros(length(in), N);
for i = 1:length(in),
    cur_in = find(comb_mat(in(i), :));
    indices(i, :) = in_from(cur_in);
end
