... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
classdef var
    % VAR represents a variable or a parameter
    
    properties (SetAccess = private, GetAccess = public)
        name = '';
        lowerBound = 0;
        upperBound = 0;
        units = 'none';
        description = '';
    end % properties
    
    properties
        value = 0;
    end % properties
    
    methods
        function v = var(name, value, lowerBound, upperBound, units, description)
            if nargin == 5
                description = '';
            end
            checkType(name, 'char', 'name');
            if isempty(regexp(name, '^[a-zA-Z_][a-zA-Z_0-9]*$'))
                error (['Name must start with a letter and be alphanumeric.  Invalid: ', name]);
            end
            if (strcmp(name, 'values') || strcmp(name, 'variables'))
                error (['Reserved name ', name]);
            end
            v.name = name;
            v.value = value;
            v.lowerBound = lowerBound;
            v.upperBound = upperBound;
            v.units = units;
            v.description = description;
        end
        
        function s = addTo(v, s)
            s{length(s) + 1} = v;
        end
        
        function s = toString(v)
            s = [v.name, ':'];
            s = pad(s, 20);
            s = sprintf('%s %.5g %s in [%.5g, %.5g]', s, v.value, v.units, v.lowerBound, v.upperBound);
            s = pad(s, 60);
            s = sprintf('%s (%s)', s, v.description);
        end
        
        function disp(v)
            fprintf(toString(v));
        end
        
    end % methods
    
end % classdef
