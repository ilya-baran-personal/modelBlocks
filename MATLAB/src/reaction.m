... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
classdef (Abstract) reaction
    %REACTION Represents a reaction in the model
    %   May be part of one or more differential equations
    properties (SetAccess = private, GetAccess = public)
        name
    end % properties

    methods
        function r = reaction(name)
            checkType(name, 'char', 'name');
            r.name = name;
        end
    end % methods

    methods (Abstract)
        reaction = prepareReaction(reaction, variables, parameters)
        modifiedDerivatives = apply(reaction, derivatives, variables, parameters)
    end % methods (Abstract)
end % classdef

