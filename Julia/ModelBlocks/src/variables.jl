# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

struct Variable
    name::String
    defaultValue::Float64
    range::Tuple{Float64, Float64}
    units::String # For now it's a string
    description::String

    function Variable(name, defaultValue, range, units, description)
        if match(r"^[a-zA-Z_][a-zA-Z0-9_]*$", name) === nothing
            error("Name must be a valid identifier: $name");
        end
        if range[1] > range[2]
            error("Invalid range in $name: $(range[1]) > $(range[2])");
        end
        new(name, defaultValue, range, units, description)
    end
end

struct Variables
    values::Vector{Float64}
    variables::Array{Variable,1}
    nameToIndex::Dict{String,Int32}
end

function Variables(variables::Array{Variable,1})
    values = [variable.defaultValue for variable = variables]
    nameToIndex = Dict{String,Int32}( variables[i].name => i for i = 1:length(variables));
    if length(nameToIndex) != length(variables)
        error("Duplicate variable names");
    end
    return Variables(values, variables, nameToIndex);
end

function Variables(variables::Variables, values::Vector{Float64})
    if length(variables.values) != length(values)
        error("Cannot create $(length(variable.values)) variables with a vector of length $(length(values))");
    end
    return Variables(values, variables.variables, variables.nameToIndex);
end

function Base.getproperty(variables::Variables, symbol::Symbol)
    if symbol === :values || symbol === :nameToIndex || symbol === :variables
        return getfield(variables, symbol)
    end
    variables.values[variables.nameToIndex[String(symbol)]]
end

function Base.getindex(variables::Variables, name::String)
    variables.values[variables.nameToIndex[name]]
end

function Base.setproperty!(variables::Variables, symbol::Symbol, value::Number)
    variables.values[variables.nameToIndex[String(symbol)]] = value;
end

function Base.setindex!(variables::Variables, value::Number, name::String)
    variables.values[variables.nameToIndex[name]] = value
end

