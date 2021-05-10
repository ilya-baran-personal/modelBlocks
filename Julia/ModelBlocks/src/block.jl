# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using OrdinaryDiffEq
using MAT

abstract type AbstractBlock end

function runBlock(block::AbstractBlock, timeRange::AbstractRange)
    floatInterval::Tuple{Float64, Float64} = convert(Tuple{Float64, Float64}, (minimum(timeRange), maximum(timeRange)));
    problem = ODEProblem((x, p, t) -> computeDerivatives(block, t, x), Vector{Float64}(getVariables(block).values), floatInterval);
    #solution = solve(problem, AutoTsit5(Rodas4()));
    solution = solve(problem, Rodas4P());
end

function solutionToVariables(solution, block::AbstractBlock, timeRange::AbstractRange)
    resultArray = Array{Array{Float64, 1}, 1}(undef, length(getVariables(block).variables));
    for index = 1:length(resultArray)
        resultArray[index] = solution(timeRange, idxs=index).u;
    end
    return Variables(getVariables(block), resultArray);
end

function solutionToMatrix(solution, timeRange::AbstractRange)::Matrix
    return Matrix(Array(solution(timeRange)));
end

function solutionToMatlab(solution, timeRange::AbstractRange, filename::String)
    t = Vector(timeRange);
    y = solutionToMatrix(solution, timeRange);
    status = string(solution.retcode);
    file = matopen(filename, "w");
    write(file, "t", t);
    write(file, "y", y);
    write(file, "status", status);
    close(file);
end

struct Block <: AbstractBlock
    variables::Variables
    parameters::Variables
    reactions::Array{AbstractReaction}
    Block(v,p,r) = new(deepcopy(v), deepcopy(p), deepcopy(r));
end

function computeDerivatives(block::Block, t::Number, x::Vector)::Vector
    derivatives = Variables(block.variables, zero(t * x));
    variables = Variables(block.variables, x);
    extraVariables = Dict{String, Any}();
    for reaction = block.reactions
        apply!(reaction, derivatives, t, variables, block.parameters, extraVariables);
    end
    #println("Time = $t");#, x = $x, result = $(derivatives.values)");
    vType = typeof(t * x);
    return vType(derivatives.values);
end

function getVariables(block::Block)::Variables block.variables; end
function getParameters(block::Block)::Variables block.parameters; end
function setParameter!(block::Block, name::String, value)
    block.parameters[name] = value;
end
function setParameters!(block::Block, values)
    block.parameters.values = values;
end

struct BlockWithBindings <: AbstractBlock
    subblock::AbstractBlock
    parameterToBinding::Dict{String, Function} # Function takes t, existing parameters, and variables and returns a value
    parameters::Variables
end

function BlockWithBindings(subblock::AbstractBlock, parameterToBinding::AbstractDict{String, T}) where T
    subblock = deepcopy(subblock);
    subblockParameters::Variables = getParameters(subblock);
    blockParameterToBinding = Dict{String, Function}();
    for (parameter, binding) in parameterToBinding
        if !haskey(subblockParameters.nameToIndex, parameter)
            error("Cannot bind nonexistent parameter $(parameter)");
        end
        if isa(binding, Function)
            blockParameterToBinding[parameter] = binding;
        else
            setParameter!(subblock, parameter, binding);
        end
    end
    newParameters = variablesSubtract(subblockParameters, Set(keys(parameterToBinding)));
    return BlockWithBindings(subblock, blockParameterToBinding, newParameters);
end

# Keep the parameters that are provided, and bind the rest to their current values
function BlockWithBindings(subblock::AbstractBlock, parametersToKeep::Array{String})
    subblock = deepcopy(subblock);
    subblockParameters::Variables = getParameters(subblock);
    newParameters::Array{Variable} = [];
    newValues::Array = [];
    for parameter in parametersToKeep
        if !haskey(subblockParameters.nameToIndex, parameter)
            error("Cannot keep nonexistent parameter $(parameter)");
        end
        index = subblockParameters.nameToIndex[parameter];
        variable = subblockParameters.variables[index];
        push!(newParameters, variable);
        push!(newValues, subblockParameters.values[index]);
    end
    newParametersV = Variables(newParameters);
    newParametersV.values = newValues;
    return BlockWithBindings(subblock, Dict{String, Function}(), newParametersV);
end

function getVariables(block::BlockWithBindings)::Variables getVariables(block.subblock); end
function getParameters(block::BlockWithBindings)::Variables block.parameters; end
function getOutputs(blockWithBindings::BlockWithBindings, timeRange::AbstractRange) getOutputs(blockWithBindings.subblock, timeRange); end
function setParameter!(block::BlockWithBindings, name::String, value)
    block.parameters[name] = value;
    setParameter!(block.subblock, name, value);
end
function setParameters!(block::BlockWithBindings, values)
    block.parameters.values = values;
    for (name, index) in block.parameters.nameToIndex
        setParameter!(block.subblock, name, block.parameters.values[index]);
    end
end

function computeDerivatives(block::BlockWithBindings, t::Number, x::Vector)::Vector
    blockVariables = Variables(getVariables(block), x);
    for (parameter, binding) in block.parameterToBinding
        setParameter!(block.subblock, parameter, binding(t, block.parameters, blockVariables));
    end
    return computeDerivatives(block.subblock, t, x);
end

struct BlockWithOutputs <: AbstractBlock
    subblock::AbstractBlock
    outputs::Variables
    computeOutputs::Function # Takes variables, parameters, solution (in Variables form), and outputs, and returns outputs
end

function getOutputs(blockWithOutputs::BlockWithOutputs, timeRange::AbstractRange)
    solution = runBlock(blockWithOutputs.subblock, timeRange);
    return blockWithOutputs.computeOutputs(getVariables(blockWithOutputs), getParameters(blockWithOutputs), timeRange,
                                           solutionToVariables(solution, blockWithOutputs, timeRange),
                                           deepcopy(blockWithOutputs.outputs));
end

function getVariables(block::BlockWithOutputs)::Variables getVariables(block.subblock); end
function getParameters(block::BlockWithOutputs)::Variables getParameters(block.subblock); end
function setParameter!(block::BlockWithOutputs, name::String, value) setParameter!(block.subblock, name, value); end
function setParameters!(block::BlockWithOutputs, values) setParameters!(block.subblock, values); end
