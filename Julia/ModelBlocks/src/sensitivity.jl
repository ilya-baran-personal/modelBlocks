# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import FiniteDiff

function localSensitivity(block::BlockWithOutputs, timeRange::AbstractRange)
    originalValues = block.block.parameters.values;
    cache = FiniteDiff.JacobianCache(originalValues, block.outputs.values);
    lambda = (x) -> begin
        block.block.parameters.values = x;
        outputs = getOutputs(block, timeRange);
        return copy(outputs.values);
    end
    jacobian = FiniteDiff.finite_difference_jacobian(lambda, originalValues);
    block.block.parameters.values = originalValues;
    return jacobian;
end
