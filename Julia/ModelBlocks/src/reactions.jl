# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

# AbstractReaction type
abstract type AbstractReaction end

function apply!(reaction::AbstractReaction, derivatives::Variables, t::Number, variables::Variables, parameters::Variables)::Nothing
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

function apply!(reaction::SimpleReaction, derivatives::Variables,  t::Number, variables::Variables, parameters::Variables)::Nothing
    rate = parameters[reaction.rateParameter];
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

# GeneralRateReaction has a function that computes the rate, then subtracts it from reactants and adds it to products
struct GeneralRateReaction <: AbstractReaction
    name::String
    reactants::Array{String}
    products::Array{String}
    rate::Function # Maps time, variables and parameters to a rate.  time is optional.
end

function apply!(reaction::GeneralRateReaction, derivatives::Variables, t::Number, variables::Variables, parameters::Variables)::Nothing
    rate = applicable(reaction.rate, variables, parameters) ? reaction.rate(variables, parameters) : reaction.rate(t, variables, parameters);
    #println("Reaction $(reaction.name) has rate $rate");
    for reactant = reaction.reactants
        derivatives[reactant] -= rate;
    end
    for product = reaction.products
        derivatives[product] += rate;
    end
end

# GeneralReaction has a function that modifies the derivatives
struct GeneralReaction <: AbstractReaction
    name::String
    apply!::Function # Modifies derivatives (first argument).  Also takes time, variables and parameters.
end

function apply!(reaction::GeneralReaction, derivatives::Variables, t::Number, variables::Variables, parameters::Variables)::Nothing
    reaction.apply!(derivatives, t, variables, parameters);
    nothing;
end

