# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import Optim
import BlackBoxOptim

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
            push!(expectedValues, outputsAndStds[variable.name][1]);
            push!(stds, outputsAndStds[variable.name][2]);
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
    return result;
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