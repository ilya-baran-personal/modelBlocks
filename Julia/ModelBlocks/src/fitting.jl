# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import Optim
import BlackBoxOptim

# Changes block parameters to the optimized values
function fitParameters!(blocksWithOutputs::AbstractArray{<:AbstractBlock},
    parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}, 1},
    outputsAndStds::AbstractDict{String, Tuple{S, T}}; # Name to value and std.
    kwargs...) where {S, T}
    result = fitParameters(blocksWithOutputs, parametersAndBounds, outputsAndStds; kwargs...);
    for block in blocksWithOutputs
        setParameters!(block, result.parameterToMinimum);
    end
    return result;
end

function fitParameters(blocksWithOutputs::AbstractArray{<:AbstractBlock},
    parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}, 1},
    outputsAndStds::AbstractDict{String, Tuple{S, T}}; # Name to value and std.
    MaxTime = 10) where {S, T}
    boundBlocks = [BlockWithBindings(block, [name for (name, lower, upper) in parametersAndBounds]) for block in blocksWithOutputs];
    lower = [lower for (_, lower, _) in parametersAndBounds];
    upper = [upper for (_, _, upper) in parametersAndBounds];

    outputs = [computeOutputs(block) for block in blocksWithOutputs];

    expectedValuesArray = [];
    stdsArray = [];
    for output in outputs
        expectedValues = [];
        stds = [];
        for variable in output.variables
            if (haskey(outputsAndStds, variable.name))
                push!(expectedValues, outputsAndStds[variable.name][1]);
                push!(stds, outputsAndStds[variable.name][2]);
            else
                push!(expectedValues, 0.0);
                push!(stds, Inf64);
            end
        end
        push!(expectedValuesArray, expectedValues);
        push!(stdsArray, stds);
    end

    objective = x -> begin
        obj = 0.;
        for i in 1:length(boundBlocks)
            setParameters!(boundBlocks[i], x);
            outputs = computeOutputs(boundBlocks[i]);
            outputs = unfold(outputs.values - expectedValuesArray[i]) ./ unfold(stdsArray[i]);
            obj += sum(outputs .^ 2)
        end
        return obj
    end

    result = BlackBoxOptim.bboptimize(objective; SearchRange = collect(zip(lower, upper)), MaxTime = MaxTime);
    parameterToMinimum = Dict{String, Any}();
    for i = 1:length(lower)
        parameterToMinimum[parametersAndBounds[i][1]] = BlackBoxOptim.best_candidate(result)[i];
    end
    return (rawResult = result, minimizer = BlackBoxOptim.best_candidate(result), parameterToMinimum = parameterToMinimum);
end

function unfold(A)
    V = []
    for x in A
        if x === A
            push!(V, x)
        else
            append!(V, unfold(x))
        end
    end
    V
end