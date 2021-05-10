# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import BlackBoxOptim

function generateVPop(blockWithOutputs::AbstractBlock,
    timeRange::AbstractRange,
    parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
    outputsAndStds::AbstractDict{String, Tuple{Float64, Float64}}, # Name to value and std
    count::Int;
    kwargs...)
    generateVPop([blockWithOutputs], timeRange, parametersAndBounds, outputsAndStds, count; kwargs...)
end

function generateVPop(blocksWithOutputs::Vector{T},
                      timeRange::AbstractRange,
                      parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
                      outputsAndStds::AbstractDict{String, Tuple{Float64, Float64}}, # Indexed to match blocks.  Name to value and std.
                      count::Int;
                      MaxTime = 1000) where T <: AbstractBlock
    boundBlocks = [BlockWithBindings(block, [name for (name, lower, upper) in parametersAndBounds]) for block in blocksWithOutputs];
    lower = [lower for (name, lower, upper) in parametersAndBounds];
    upper = [upper for (name, lower, upper) in parametersAndBounds];

    outputs = [getOutputs(block, timeRange) for block in blocksWithOutputs];

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
            outputs = getOutputs(boundBlocks[i], timeRange);
            outputs = (outputs.values - expectedValuesArray[i]) ./ stdsArray[i];
            obj += sum(max.(outputs .^ 2, 1) .- 1)
        end
        return obj
    end

    vpop = Matrix{Float64}(undef, length(lower), count);

    for i in 1:count
        result = BlackBoxOptim.bboptimize(objective; SearchRange = collect(zip(lower, upper)), MaxTime = MaxTime / count);
        vpop[:, i] = BlackBoxOptim.best_candidate(result);
    end
    vpop;
end