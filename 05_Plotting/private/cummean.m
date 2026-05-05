function [mu,n] = cummean(data, group)
    tmp = splitapply(@(x){cumsum(x(:), "omitmissing")}, data, group);
    

    n = splitapply(@(x){cumsum(~isnan(x(:)))}, data, group);
    n = [n{:}].';

    mu = [tmp{:}].' ./ n;
end