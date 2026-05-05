function [sigma, n] = cumstd(data, group)
    EX2 = splitapply(@(x){cumsum(x(:).^2, "omitmissing")}, data, group);
    EX2 = [EX2{:}].';

    n = splitapply(@(x){cumsum(~isnan(x(:)))}, data, group);
    n = [n{:}].';

    mu = cummean(data, group);
    
    sigma = sqrt( (EX2 - n .* mu.^2) ./ (n - 1) );
end