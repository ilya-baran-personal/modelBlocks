# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import BlackBoxOptim

function generateVPop(blockWithOutputs::AbstractBlock,
                      timeRange::AbstractRange,
                      parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}},
                      outputsAndStds::AbstractDict{String, Tuple{Float64, Float64}}, # Name to value and std
                      count::Int;
                      MaxTime = 1000)
    boundBlock = BlockWithBindings(blockWithOutputs, [name for (name, lower, upper) in parametersAndBounds]);
    lower = [lower for (name, lower, upper) in parametersAndBounds];
    upper = [upper for (name, lower, upper) in parametersAndBounds];

    outputs = getOutputs(blockWithOutputs, timeRange);

    expectedValues = []
    stds = [];
    for variable in outputs.variables
        push!(expectedValues, outputsAndStds[variable.name][1]);
        push!(stds, outputsAndStds[variable.name][2]);
    end

    objective = x -> begin    
        setParameters!(blockWithOutputs, x);
        outputs = getOutputs(blockWithOutputs, timeRange);
        outputs = (outputs.values - expectedValues) ./ stds;
        obj = sum(max.(outputs .^ 2, 1) .- 1)
        return obj
    end

    vpop = Matrix{Float64}(undef, length(getParameters(blockWithOutputs).variables), count);

    for i in 1:count
        result = BlackBoxOptim.bboptimize(objective; SearchRange = collect(zip(lower, upper)), MaxTime = MaxTime / count);
        vpop[:, i] = BlackBoxOptim.best_candidate(result);
    end
    vpop;
end