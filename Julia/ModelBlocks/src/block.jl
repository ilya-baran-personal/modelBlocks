# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using OrdinaryDiffEq
using MAT

abstract type AbstractBlock end

struct BlockOutputDefinition
    outputs::Variables
    computeOutputs::Function # Takes variables, parameters, solution (in Variables form), and outputs, and returns outputs
end

mutable struct BlockExtraData
    timeRange::Union{AbstractRange, Nothing}
    discontinuities::Vector
    outputDefinition::Union{BlockOutputDefinition, Nothing}
    BlockExtraData() = new(nothing, [], nothing)
end

struct Block <: AbstractBlock
    variables::Variables
    parameters::Variables
    reactions::Array{AbstractReaction}
    extraData::BlockExtraData
    Block(v, p, r) = new(deepcopy(v), deepcopy(p), deepcopy(r), BlockExtraData());
end

function getExtraData(block::AbstractBlock)::BlockExtraData block.extraData end

function getTimeRange(block::AbstractBlock)::AbstractRange
    range = getExtraData(block).timeRange;
    if range === nothing
        error("Expected time range but it was missing");
    end
    range;
end

function setTimeRange!(block::AbstractBlock, timeRange::AbstractRange) getExtraData(block).timeRange = timeRange; end

function getOutputDefinition(block::AbstractBlock)::BlockOutputDefinition
    outputs = getExtraData(block).outputDefinition;
    if outputs === nothing
        error("Expected output definition but it was missing");
    end
    outputs;
end

function setOutputDefinition!(block::AbstractBlock, outputs::Variables, computeOutputs::Function)
    getExtraData(block).outputDefinition = BlockOutputDefinition(outputs, computeOutputs);
end

function getDiscontinuities(block::AbstractBlock)::Vector getExtraData(block).discontinuities; end

function setDiscontinuities!(block::AbstractBlock, discontinuities::Vector{T}) where T <: Number
    getExtraData(block).discontinuities = discontinuities;
end

function runBlock(block::AbstractBlock)
    timeRange = getTimeRange(block);
    floatInterval::Tuple{Float64, Float64} = convert(Tuple{Float64, Float64}, (minimum(timeRange), maximum(timeRange)));
    problem = ODEProblem((x, p, t) -> computeDerivatives(block, t, x), Vector{Float64}(getVariables(block).values), floatInterval);
    solve(problem, Rodas4P(), tstops = getDiscontinuities(block));
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

function computeOutputs(blockWithOutputs::AbstractBlock)
    outputDefinition = getOutputDefinition(blockWithOutputs);
    timeRange = getTimeRange(blockWithOutputs);
    solution = runBlock(blockWithOutputs);
    return outputDefinition.computeOutputs(getVariables(blockWithOutputs), getParameters(blockWithOutputs), timeRange,
                                           solutionToVariables(solution, blockWithOutputs, timeRange),
                                           deepcopy(outputDefinition.outputs));
end

# ============================================ BlockWithBindings =======================================
struct BlockWithBindings <: AbstractBlock
    subblock::AbstractBlock
    parameterToBinding::Dict{String, Function} # Function takes t, variables, and existing parameters and returns a value
    parameters::Variables
    extraData::BlockExtraData
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
    extraData = BlockExtraData();
    extraData.timeRange = getExtraData(subblock).timeRange;
    return BlockWithBindings(subblock, blockParameterToBinding, newParameters, extraData);
end

# Keep the parameters that are provided, and bind the rest to their current values
function BlockWithBindings(subblock::AbstractBlock, parametersToKeep::AbstractArray{String})
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
    extraData = BlockExtraData();
    extraData.timeRange = getExtraData(subblock).timeRange;
    return BlockWithBindings(subblock, Dict{String, Function}(), newParametersV, extraData);
end

function getVariables(block::BlockWithBindings)::Variables getVariables(block.subblock); end
function getParameters(block::BlockWithBindings)::Variables block.parameters; end
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
        setParameter!(block.subblock, parameter, binding(t, blockVariables, block.parameters));
    end
    return computeDerivatives(block.subblock, t, x);
end

function computeOutputs(blockWithOutputs::BlockWithBindings)
    timeRange = getTimeRange(blockWithOutputs);
    solution = runBlock(blockWithOutputs);

    outputDefinition = getExtraData(blockWithOutputs).outputDefinition;
    if (outputDefinition === nothing)
        outputDefinition = getOutputDefinition(blockWithOutputs.subblock);
        return outputDefinition.computeOutputs(getVariables(blockWithOutputs), getParameters(blockWithOutputs.subblock), timeRange,
                                               solutionToVariables(solution, blockWithOutputs, timeRange),
                                               deepcopy(outputDefinition.outputs));
    end
    return outputDefinition.computeOutputs(getVariables(blockWithOutputs), getParameters(blockWithOutputs), timeRange,
                                           solutionToVariables(solution, blockWithOutputs, timeRange),
                                           deepcopy(outputDefinition.outputs));
end

# ============================================ BlockCombo =======================================

const EXTRA_NAME = "_extra_";

struct BlockCombo <: AbstractBlock
    subblocks::Vector{Tuple{String, AbstractBlock}}
    parameters::Variables
    parameterNameToBlocksAndNames::Dict{String, Vector{Tuple{String, String}}} # For setting parameters
    glueFunctions::Vector{Tuple{String, String, Function}} # Block name, parameter name, function that computes it
    # Glue functions take time, variables, parameters (of the entire block)
    variables::Variables
    extraData::BlockExtraData
end

function BlockCombo(subblocks::Vector{Tuple{String, T}},
                    glueFunctions::Vector{Tuple{String, String, F}},
                    extraParameters::Variables) where {T <: AbstractBlock, F <: Function}
    if (length(subblocks) < 1)
        error("Must have at least one subblock");
    end
    extraData = getExtraData(subblocks[1][2]);
    subblocks = deepcopy(subblocks);
    gluedParameters = Set((f[1], f[2]) for f in glueFunctions);
    
    # Gather parameters
    seenNames = Set{String}();
    parameterArray = Vector{Variable}();
    parameterValueArray = [];
    parameterNameToBlocksAndNames = Dict{String, Vector{Tuple{String, String}}}();

    namesAndParameterSets = [(b[1], getParameters(b[2])) for b in subblocks];
    push!(namesAndParameterSets, (EXTRA_NAME, extraParameters));

    for namesAndParameters in namesAndParameterSets
        parameters = namesAndParameters[2];
        for parameter::Variable in parameters.variables
            if (in((namesAndParameters[1], parameter.name), gluedParameters))
                continue;
            end
            if !in(parameter.name, seenNames)
                push!(seenNames, parameter.name);
                push!(parameterArray, parameter);
                push!(parameterValueArray, parameters.values[parameters.nameToIndex[parameter.name]]);                
            end
            if (namesAndParameters[1] != EXTRA_NAME)
                push!(get!(parameterNameToBlocksAndNames, parameter.name, []), (namesAndParameters[1], parameter.name));
            end
        end
    end

    parameters = Variables(parameterArray);

    variableArray = Vector{Variable}();
    variableValueArray = [];

    for nameAndSubblock in subblocks
        variables = getVariables(nameAndSubblock[2]);
        for variable in variables.variables
            push!(variableArray, variable);
            push!(variableValueArray, variables.values[variables.nameToIndex[variable.name]]);
        end
    end

    variables = Variables(variableArray);
    variables.values = variableValueArray;

    result = BlockCombo(subblocks, parameters, parameterNameToBlocksAndNames, glueFunctions, variables, extraData);
    setParameters!(result, parameterValueArray); # Propagate to subblocks
    result;
end

function getSubblock(block::BlockCombo, name::String)
    for (subblockName, subblock) in block.subblocks
        if name == subblockName
            return subblock;
        end
    end
    error("Subblock $name not found");
end

function getVariables(block::BlockCombo)::Variables block.variables; end
function getParameters(block::BlockCombo)::Variables block.parameters; end
function setParameter!(block::BlockCombo, name::String, value)
    block.parameters[name] = value;
    blocksAndNames = block.parameterNameToBlocksAndNames[name];
    for blockAndName in blocksAndNames
        setParameter!(getSubblock(block, blockAndName[1]), blockAndName[2], value);
    end
end

function setParameters!(block::BlockCombo, values)
    block.parameters.values = values;
    
    for parameter in block.parameters.variables
        blocksAndNames = get(block.parameterNameToBlocksAndNames, parameter.name, nothing);
        if blocksAndNames === nothing
            continue;
        end
        for blockAndName in blocksAndNames
            setParameter!(getSubblock(block, blockAndName[1]), blockAndName[2], values[block.parameters.nameToIndex[parameter.name]]);
        end
    end
end

function computeDerivatives(block::BlockCombo, t::Number, x::Vector)::Vector
    blockVariables = Variables(block.variables, x);
    for nameParameterAndGlue in block.glueFunctions
        value = nameParameterAndGlue[3](t, blockVariables, block.parameters);
        setParameter!(getSubblock(block, nameParameterAndGlue[1]), nameParameterAndGlue[2], value);
    end

    derivatives = [];
    variableIndex = 1;
    for (_, subblock) in block.subblocks
        nVariables = length(getVariables(subblock).variables);
        push!(derivatives, computeDerivatives(subblock, t, x[variableIndex:(variableIndex + nVariables - 1)]))
        variableIndex += nVariables;
    end

    vcat(derivatives...);
end
