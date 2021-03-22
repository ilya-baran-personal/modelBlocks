# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

# AbstractReaction type
abstract type AbstractReaction end

function apply!(reaction::AbstractReaction, derivatives::Variables, variables::Variables, parameters::Variables)::Nothing
    error("Function apply! must be defined for concrete reactions");
end

function reactionName(reaction::AbstractReaction)::String reaction.name end

# SimpleReaction is of the form aA + bB + ... --> cC + dD + ...
struct SimpleReaction <: AbstractReaction
    name::String
    reactants::Array{String}
    products::Array{String}
    rateParameter::String
end

function apply!(reaction::SimpleReaction, derivatives::Variables, variables::Variables, parameters::Variables)::Nothing
    rate::Float64 = parameters[reaction.rateParameter];
    for reactant = reaction.reactants
        rate *= variables[reactant];
    end
    #println("Reaction $(reaction.name) has rate $rate and parameter $(parameters[reaction.rateParameter])");
    for reactant = reaction.reactants
        derivatives[reactant] -= rate;
    end
    for product = reaction.products
        derivatives[product] += rate;
    end
end
