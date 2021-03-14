... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
classdef block
    %BLOCK Represents a model block
    %   Has parameters, variables, reactions, dxdt function
    
    properties
        variables
        parameters
        reactions        
    end % properties
    
    methods
        function b = block(v, p)
            b.variables = v;
            b.parameters = p;
        end
    end % methods
end % classdef
