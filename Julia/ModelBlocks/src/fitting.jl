# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import Optim
import BlackBoxOptim

function fitParameters(block::BlockWithOutputs, timeRange::AbstractRange, expectedOutputs::Variables)
    lower = [variable.range[1] for variable in getParameters(block.block).variables];
    upper = [variable.range[2] for variable in getParameters(block.block).variables];
    initial_x = getParameters(block.block).values;
    lower = min.(initial_x * 0.5, initial_x * 2);
    upper = max.(initial_x * 0.5, initial_x * 2);
    return fitParameters(block, timeRange, lower, upper, expectedOutputs)
end

function fitParameters(block::BlockWithOutputs, timeRange::AbstractRange,
                       parametersAndBounds::AbstractDict{String, Tuple{Float64, Float64}}, expectedOutputs::Variables)
    boundBlock = BlockWithBindings(block.block, collect(keys(parametersAndBounds)));
    lower = [value[1] for (key, value) in parametersAndBounds];
    upper = [value[2] for (key, value) in parametersAndBounds];
    boundWithOutputs = BlockWithOutputs(boundBlock, block.outputs, block.computeOutputs);
    return fitParameters(boundWithOutputs, timeRange, lower, upper, expectedOutputs);
end

function fitParameters(block::BlockWithOutputs, timeRange::AbstractRange, lower::Array, upper::Array, expectedOutputs::Variables)
    initial_x = deepcopy(getParameters(block.block).values);
    count = 1;
    f = x -> begin
        setParameters!(block.block, x);
        outputs = getOutputs(block, timeRange);
        difference = outputs.values - expectedOutputs.values;
        err = sum(difference .^ 2);
        if mod(count += 1, 100) == 0
            println("c = $count, x = $x, err = $err");
        end
        return err;
    end

    #result = Optim.optimize(f, lower, upper, initial_x, SAMIN(rt=0.5), Optim.Options(iterations=10^4));
    result = BlackBoxOptim.bboptimize(f; SearchRange = collect(zip(lower, upper)), MaxTime = 10.0);
    setParameters!(block.block, initial_x);
    return result;
end
