... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
classdef block
    %BLOCK Represents a model block
    %   Has parameters, variables, reactions, dxdt function
    
    properties
        variables
        parameters
        reactions = {}
    end % properties
    
    methods
        function b = block(variables, parameters)
            checkType(variables, 'vars', 'variables');
            checkType(parameters, 'vars', 'parameters');
            b.variables = variables;
            b.parameters = parameters;
        end
        
        function b = addReaction(b, r)
            checkType(r, 'reaction', 'input to addReaction');
            b.reactions{numel(b.reactions) + 1} = r;
        end
        
        function derivativesVector = computeDerivatives(b, t, variablesVector, p)
            derivativesVector = zeros(size(variablesVector));
            for i = 1:numel(b.reactions)
                derivativesVector = b.reactions{i}.apply(derivativesVector, variablesVector, p);
            end
        end
        
        function odefun = odefun(b, p)
            odefun = @(t, v) b.computeDerivatives(t, v, p);
        end
        
        function solve = run(b, tInterval)
            for i = 1:numel(b.reactions)
                b.reactions{i} = b.reactions{i}.prepareReaction(b.variables, b.parameters);
            end
            solve = ode15s(b.odefun(b.parameters.values), tInterval, b.variables.values);
        end
    end % methods
end % classdef
