function [train_err, test_err] = calculate_error (D, train_patterns, train_targets, test_patterns, test_targets, region, Nclasses)

% Calculate error (used by the main calculation functions)

train_err   = zeros(Nclasses+1,1);
test_err    = zeros(Nclasses+1,1);

if ~isempty(train_targets),    
    [classify, err]  = classification_error(D, train_patterns, train_targets, region);
    for j = 1:Nclasses,
        train_err(j)   = 1 - classify(j,j);    
    end
    if (Nclasses>=2)
        train_err(Nclasses+1) = err;     
    else
        train_err(3) = err;
    end
end

if ~isempty(test_targets),
    [classify, err]  = classification_error(D, test_patterns, test_targets, region);
    for j = 1:Nclasses,
       test_err(j)   = 1 - classify(j,j);    
    end
    if (Nclasses>=2)
        test_err(Nclasses+1) = err;     
    else
        test_err(3) = err;
    end
end
