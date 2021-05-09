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

Base.show(io::IO, v::Variable) = print(io, "Var $(v.name) in $(v.range) default: $(v.defaultValue) units: $(v.units) $(v.description)");

struct Variables
    values::Vector{Any}
    variables::Array{Variable,1}
    nameToIndex::Dict{String,Int32}
end

function Variables(variables::Array{Variable,1})
    values = Vector{Any}([variable.defaultValue for variable = variables])
    nameToIndex = Dict{String,Int32}(variables[i].name => i for i = 1:length(variables));
    if length(nameToIndex) != length(variables)
        error("Duplicate variable names");
    end
    return Variables(values, copy(variables), nameToIndex);
end

function Variables(variables::Variables, values::Vector{Any})
    if length(variables.values) != length(values)
        error("Cannot create $(length(variables.values)) variables with a vector of length $(length(values))");
    end
    return Variables(values, variables.variables, variables.nameToIndex);
end

function Variables(variables::Variables, values::Vector)
    return Variables(variables, Vector{Any}(values));
end

function variablesToMatlab(variables::Variables, filename::String)
    v = Dict{String,Any}(variables.variables[i].name => variables.values[i] for i = 1:length(variables.values));
    file = matopen(filename, "w");
    write(file, "v", v);
    close(file);
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

function Base.setproperty!(variables::Variables, symbol::Symbol, value)
    if symbol === :values
        variables.values[1:end] = value;
        return
    end
    variables.values[variables.nameToIndex[String(symbol)]] = value;
end

function Base.setindex!(variables::Variables, value, name::String)
    variables.values[variables.nameToIndex[name]] = value
end

function Base.show(io::IO, ::MIME"text/plain", vars::Variables)
    println("Variables:")
    maxLen = maximum(v -> length(v.name), vars.variables);
    for i = 1:length(vars.variables)
        v = vars.variables[i];
        value = vars.values[i];
        println(io, "$(lpad(v.name, maxLen + 2)) = $value $(v.units) in $(v.range) $(v.description)");
    end
end

# Combining variables
function variablesUnion(v1::Variables, v2::Variables)::Variables
    newSize = length(v1.variables);
    mergedNameToIndex = copy(v1.nameToIndex);
    for var in v2.variables
        if !haskey(mergedNameToIndex, var.name)
            mergedNameToIndex[var.name] = newSize += 1;
        end
    end
    newVariables = Array{Variable}(undef, newSize);
    newValues = Vector{Any}(undef, newSize);
    newVariables[1:length(v1.variables)] = v1.variables;
    newValues[1:length(v1.variables)] = v1.values;
    for (name, index) in v2.nameToIndex
        newIndex = mergedNameToIndex[name];
        newValues[newIndex] = v2.values[index];
        newVariables[newIndex] = v2.variables[index];
    end
    return Variables(newValues, newVariables, mergedNameToIndex);
end

function variablesSubtract(v::Variables, toRemove::AbstractSet{String})::Variables
    newValues = Vector{Any}();
    newVariables = Array{Variable,1}();
    newNameToIndex = Dict{String,Int32}();    
    for var in v.variables
        if in(var.name, toRemove)
            continue;
        end
        index = v.nameToIndex[var.name];
        push!(newVariables, var);
        push!(newValues, v.values[index]);
        newNameToIndex[var.name] = length(newValues);
    end
    return Variables(newValues, newVariables, newNameToIndex);
end