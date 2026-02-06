classdef DataLogger < handle
    % Is responsible for caching and saving user-relevant data
   
    properties (Access = private)
        file_path   (1,1) string
        cap         (1,1) double
        data        (1,1) struct    = struct() % fields: varname -> cell buffer
        count       (1,1) struct    = struct() % fields: varname -> number of logged entries
        linked_log                  = @(varargin) []
        linked_delete               = @(varargin) []
    end

    methods
        function obj = DataLogger(params)
            arguments
                params.log_enabled      (1,1) logical               = false
                params.add_timestamp    (1,1) logical               = false
                params.file_location    (1,1) string {mustBeFolder} = "02_Data"
                params.file_name        (1,1) string                = ""
                params.initial_capacity (1,1) double                = 1000
            end
            obj.cap         = params.initial_capacity;
            if params.log_enabled
                obj.linked_log      = @obj.log_implementation;
                obj.linked_delete   = @obj.delete_implementation;

                if params.file_name == ""
                    error("DataLogger:FileName", "Filename cannot be empty");
                end
            end

            obj.file_path   = params.file_location + "\" + params.file_name;
            if params.add_timestamp
                obj.file_path = obj.file_path + "_" + string(datetime('now', 'Format','yyyy-MM-dd_HHmmss'));
            end
            

        end

        function log(obj, varargin)
            obj.linked_log(varargin{:});
        end

        function delete(obj)
            obj.linked_delete();
        end

    end

    methods (Access = private)

        function log_implementation(obj, varargin)
            if mod(nargin-1, 2) ~= 0
                error("DataLogger:BadArgs", "Use name-value pairs: log(name1, val1, name2, val2, ...)");
            end

            for i = 1:2:numel(varargin)
                name    = obj.normalize_name(varargin{i});
                val     = varargin{i+1};

                if ~isfield(obj.data, name)
                    obj.data.(name) = cell(obj.cap, 1);
                    obj.count.(name) = 0;
                end

                obj.count.(name) = obj.count.(name) + 1;
                k = obj.count.(name);

                if k > numel(obj.data.(name))
                    obj.grow_field(name);
                end

                obj.data.(name){k,1} = val;
            end
        end
        
        function delete_implementation(obj)
            try
                S = struct();
                f = fieldnames(obj.data);

                for idx = 1:numel(f)
                    name = f{idx};
                    k  = obj.count.(name);
                    S.(name) = obj.data.(name)(1:k,1);
                    S.(name) = [S.(name){:}];
                end

                save(obj.file_path, "-struct", "S");
                fprintf("Saved data in %s.mat\n", obj.file_path);
            catch ME
                warning("DataLogger:SaveFailed", "Could not save log: %s", ME.message);
            end
        end

        function grow_field(obj, name)
            old_buf = obj.data.(name);
            new_cap = max( 2*numel(old_buf), numel(old_buf) + obj.cap );
            old_buf(new_cap,1) = {[]}; 
            obj.data.(name) = old_buf;
        end

        function name = normalize_name(~, name)
            if isstring(name), name = char(name); end
            if ~(ischar(name) && ~isempty(name))
                error("DataLogger:BadName", "Names must be char or string scalars.");
            end
            name = matlab.lang.makeValidName(name);
        end
    end
end