... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
classdef basicReaction < reaction
    %BASICREACTION A reaction of the form aA+bB+... -> cC+dD+...
    %   Rate is proportional to the product of reactant concentrations (raised to
    %   stoichiometric coefficient powers)
    
    properties (SetAccess = private, GetAccess = public)
        reactants % cell array of strings
        products % cell array of strings
        rateParameter % string
        reactantIndices;
        productIndices;
        rateParameterIndex;
    end
    
    methods
        function obj = basicReaction(name, reactants, products, rateParameter)
            obj = obj@reaction(name);
            checkType(reactants, 'cell', 'reactants');
            obj.reactants = reactants;
            checkType(products, 'cell', 'products');
            obj.products = products;
            checkType(rateParameter, 'char', 'rateParameter');
            obj.rateParameter = rateParameter;
        end
        
        function reaction = prepareReaction(reaction, variables, parameters)
            reaction.reactantIndices = zeros(numel(reaction.reactants), 1);
            for i = 1:numel(reaction.reactants)
                reaction.reactantIndices(i) = variables.(reaction.reactants{i}).index;
            end            
            reaction.productIndices = zeros(numel(reaction.products), 1);
            for i = 1:numel(reaction.products)
                reaction.productIndices(i) = variables.(reaction.products{i}).index;
            end
            reaction.rateParameterIndex = parameters.(reaction.rateParameter).index;
        end
        
        function derivatives = apply(reaction, derivatives, variables, parameters)
            % Compute rate
            rate = parameters(reaction.rateParameterIndex) * prod(variables(reaction.reactantIndices));
            % Apply changes
            derivatives(reaction.reactantIndices) = derivatives(reaction.reactantIndices) - rate;
            derivatives(reaction.productIndices) = derivatives(reaction.productIndices) - rate;
        end
    end
end

