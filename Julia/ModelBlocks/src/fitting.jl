# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using Optim

function fitParameters(block::BlockWithOutputs, timeRange::AbstractRange, expectedOutputs::Variables)
    lower = [variable.range[1] for variable in block.block.parameters.variables];
    upper = [variable.range[2] for variable in block.block.parameters.variables];
    initial_x = block.block.parameters.values;
    lower = min.(initial_x * 0.5, initial_x * 2);
    upper = max.(initial_x * 0.5, initial_x * 2);
    return fitParameters(block, timeRange, lower, upper, expectedOutputs)
end

function fitParameters(block::BlockWithOutputs, timeRange::AbstractRange, lower::Array, upper::Array, expectedOutputs::Variables)
    initial_x = block.block.parameters.values;
    count = 1;
    f = x -> begin
        block.block.parameters.values = x;
        outputs = getOutputs(block, timeRange);
        difference = outputs.values - expectedOutputs.values;
        err = sum(difference .^ 2);
        if mod(count += 1, 100) == 0
            println("c = $count, x = $x, err = $err");
        end
        return err;
    end

    result = optimize(f, lower, upper, initial_x, SAMIN(rt=0.5), Optim.Options(iterations=10^4));
    block.block.parameters.values = initial_x;
    return result;
end