# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

struct Variable
    name::String
    defaultValue::Float64
    range::Tuple{Float64, Float64}
    units::String # For now it's a string
    description::String
end

struct Variables
    values::Array{Float64,1}
    variables::Array{Variable,1}
    nameToIndex::Dict{String,Int32}
end

function Variables(variables::Array{Variable,1})
    values = [variable.defaultValue for variable = variables]
    nameToIndex = Dict{String,Int32}( variables[i].name => i for i = 1:length(variables));
    return Variables(values, variables, nameToIndex)
end

function Base.getproperty(variables::Variables, symbol::Symbol)
    if symbol === :values || symbol === :nameToIndex
        return getfield(variables, symbol)
    end
    variables.values[variables.nameToIndex[String(symbol)]]
end

function Base.setproperty!(variables::Variables, symbol::Symbol, value::Number)
    variables.values[variables.nameToIndex[String(symbol)]] = value;
end
