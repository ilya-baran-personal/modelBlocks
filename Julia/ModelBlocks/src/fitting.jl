# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using Optim

function fitParameters(block::BlockWithOutputs, timeRange::AbstractRange, expectedOutputs::Variables)
    lower = [variable.range[1] for variable in getParameters(block.block).variables];
    upper = [variable.range[2] for variable in getParameters(block.block).variables];
    initial_x = getParameters(block.block).values;
    lower = min.(initial_x * 0.5, initial_x * 2);
    upper = max.(initial_x * 0.5, initial_x * 2);
    return fitParameters(block, timeRange, lower, upper, expectedOutputs)
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

    result = optimize(f, lower, upper, initial_x, SAMIN(rt=0.5), Optim.Options(iterations=10^4));
    setParameters!(block.block, initial_x);
    return result;
end