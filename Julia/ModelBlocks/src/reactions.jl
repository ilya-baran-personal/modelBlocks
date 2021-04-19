# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

# AbstractReaction type
abstract type AbstractReaction end

function apply!(reaction::AbstractReaction, derivatives::Variables, t::Number, variables::Variables, parameters::Variables, extraVariables::Dict{String, Any})::Nothing
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

function apply!(reaction::SimpleReaction, derivatives::Variables,  t::Number, variables::Variables, parameters::Variables, extraVariables::Dict{String, Any})::Nothing
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
    rate::Function # Maps time, variables, parameters, and extraVariables to a rate.  time and extraVariables are optional, but if extraVariables are provided, time must be also
    usesT::Bool
    usesTAndExtraVariables::Bool
    function GeneralRateReaction(name, reactants, products, rate)
        variables = Variables(Array{Variable, 1}());
        usesT = !applicable(rate, variables, variables);
        usesTAndExtraVariables = usesT && !applicable(rate, 0, variables, variables);
        new(name, reactants, products, rate, usesT, usesTAndExtraVariables);
    end
end

function apply!(reaction::GeneralRateReaction, derivatives::Variables, t::Number, variables::Variables, parameters::Variables, extraVariables::Dict{String, Any})::Nothing
    if reaction.usesTAndExtraVariables
        rate = reaction.rate(t, variables, parameters, extraVariables);
    elseif reaction.usesT
        rate = reaction.rate(t, variables, parameters);
    else
        rate = reaction.rate(variables, parameters);
    end
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
    apply!::Function # Modifies derivatives (first argument).  Also takes time, variables, parameters, and extraVariables (which it can modify).
end

function apply!(reaction::GeneralReaction, derivatives::Variables, t::Number, variables::Variables, parameters::Variables, extraVariables::Dict{String, Any})::Nothing
    reaction.apply!(derivatives, t, variables, parameters, extraVariables);
    nothing;
end
