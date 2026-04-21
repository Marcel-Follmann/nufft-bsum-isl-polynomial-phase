function data_out = reorder_data(data_in, order, varargin)
    size_vec = zeros(1,nargin-2);
    for i = 1:nargin-2
        len = numel(varargin{i});
        if len == 1
            len = varargin{i};
        end
        size_vec(i) = len;
    end

    data_reshaped = reshape(data_in, size_vec);
    data_out = permute(data_reshaped, order);
end