# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import Optim
import BlackBoxOptim

function fitParameters(block::AbstractBlock, timeRange::AbstractRange, expectedOutputs::Variables; kwargs...)
    lower = [variable.range[1] for variable in getParameters(block).variables];
    upper = [variable.range[2] for variable in getParameters(block).variables];
    initial_x = getParameters(block).values;
    lower = min.(initial_x * 0.5, initial_x * 2);
    upper = max.(initial_x * 0.5, initial_x * 2);
    return fitParameters(block, timeRange, lower, upper, expectedOutputs; kwargs...)
end

function fitParameters(block::AbstractBlock, timeRange::AbstractRange,
                       parametersAndBounds::AbstractArray{Tuple{String, Float64, Float64}}, expectedOutputs::Variables; kwargs...)
    boundBlock = BlockWithBindings(block, [name for (name, lower, upper) in parametersAndBounds]);
    lower = [lower for (name, lower, upper) in parametersAndBounds];
    upper = [upper for (name, lower, upper) in parametersAndBounds];
    return fitParameters(boundBlock, timeRange, lower, upper, expectedOutputs; kwargs...);
end

function fitParameters(block::AbstractBlock, timeRange::AbstractRange, lower::Array, upper::Array, expectedOutputs::Variables; MaxTime = 10)
    initial_x = deepcopy(getParameters(block).values);
    count = 1;
    f = x -> begin
        setParameters!(block, x);
        outputs = getOutputs(block, timeRange);
        difference = unfold(outputs.values) - unfold(expectedOutputs.values);
        err = sum(difference .^ 2);
        if mod(count += 1, 100) == 0
            println("c = $count, x = $x, err = $err, diff = $difference");
        end
        return err;
    end

    #result = Optim.optimize(f, lower, upper, initial_x, SAMIN(rt=0.5), Optim.Options(iterations=10^4));
    result = BlackBoxOptim.bboptimize(f; SearchRange = collect(zip(lower, upper)), MaxTime = MaxTime);
    setParameters!(block, initial_x);    
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