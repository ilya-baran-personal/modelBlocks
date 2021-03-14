... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
classdef vars
    % VARS represents a set of variables or parameters
    % The idea is to be able to refer to variables by either name or index
    % See subsref documentation
    
    properties (Access = private)
        variables = {}
        variableToIndex
        n = 0
    end % properties
    
    methods
        function v = vars(varCollection)
            if (isstruct(varCollection))
                varCollection = orderfields(varCollection);
                names = fieldnames(varCollection);
                for i = 1:size(names, 1)
                    v = addVariable(v, varCollection.(key));
                end
            elseif (iscell(varCollection))
                for i = 1:numel(varCollection)
                    v = addVariable(v, varCollection{i});
                end
            end
        end
        
        function v = addVariable(v, var)
            if (~isa(var, 'var'))
                error ('Cannot call addVariable with something that is not a var');
            end
            if (isfield(v.variableToIndex, var.name))
                error (['Duplicate variable name ', var.name]);
            end
            v.n = v.n + 1;
            v.variableToIndex.(var.name) = v.n;
            v.variables{v.n} = var;
        end
        
        function s = toString(v)
            s = '';
            for variable = v.variables
                s = [s, variable{1}.toString(), '\n'];
            end
        end
        
        function disp(v)
            fprintf(toString(v));
        end
        
        function n = size(v)
            n = v.n;
        end
        
        % Valid indices are:
        %   .varName -- gets the value of the variable varName
        %   .varName.upperBound (and other variable properties)
        %   .varName.index -- gets the index
        %   .values -- gets a vector of variable values
        %   .variables -- gets a cell array of the variables
        function r = subsref(v, idx)
            if (isempty (idx))
                error ('@vars/subsref: missing index');
            end
            
            switch (idx(1).type)
                case '()'
                    error ('@vars/subsref: invalid index');
                case '{}'
                    error ('@vars/subsref: invalid index');
                case '.'
                    field = idx(1).subs;
                    if (isfield(v.variableToIndex, field))
                        r = v.variables{v.variableToIndex.(field)};
                        if (numel(idx) == 1)
                            r = r.value;
                        elseif (strcmp(idx(2).type, '.') && strcmp(idx(2).subs, 'index'))
                            r = v.variableToIndex.(field);
                            idx = idx(2:end);
                        end
                    else
                        switch (field)
                            case 'values'
                                r = zeros(1, v.n);
                                for i = 1:v.n
                                    r(1, i) = v.variables{i}.value;
                                end
                            case 'variables'
                                r = v.variables;
                            otherwise
                                error (['@vars/subsref: invalid field ', field]);
                        end
                    end
            end
            
            if (numel (idx) > 1)
                r = subsref (r, idx(2:end));
            end
        end % function subsref
        
        function v = subsasgn(v, idx, rhs)
            if (isempty (idx))
                error ('@vars/subsasgn: missing index');
            end
            
            switch (idx(1).type)
                case '()'
                    error ('@vars/subsasgn: invalid index');
                case '{}'
                    error ('@vars/subsasgn: invalid index');
                case '.'
                    field = idx(1).subs;
                    if (isfield(v.variableToIndex, field))
                        if (numel(idx) > 1)
                            error ('@vars/subsasgn: too many indices');
                        end
                        
                        v.variables{v.variableToIndex.(field)}.value = rhs;
                    elseif (strcmp('values', field))
                        if (numel(idx) == 1) % Assign the whole values vector
                            if (numel(rhs) ~= v.n)
                                error ('Expected %d values, got %d', v.n, numel(rhs));
                            end
                            values = reshape(rhs, 1, v.n);
                        else % Allow subindexing into values
                            values = zeros(1, v.n);
                            for i = 1:v.n
                                values(1, i) = v.variables{i}.value;
                            end
                            values = subsasgn(values, idx(2:end), rhs);
                        end
                        for i = 1:v.n
                            v.variables{i}.value = values(1, i);
                        end
                    else
                        error ('@vars/subsasgn: invalid index');
                    end
            end % switch
        end % function subsasgn
    end % methods
    
    methods (Hidden)
        function checkProps(v) % for debugging private properties
            keyboard
        end
    end % methods (Hidden)
end % classdef
