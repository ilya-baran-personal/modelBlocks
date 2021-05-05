# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import FiniteDiff

function localSensitivity(block::BlockWithOutputs, timeRange::AbstractRange)
    originalValues = Vector{Float64}(getParameters(block.block).values);
    lambda = (x) -> begin
        setParameters!(block.block, x);
        outputs = getOutputs(block, timeRange);
        return Vector{Float64}(outputs.values);
    end
    jacobian = FiniteDiff.finite_difference_jacobian(lambda, originalValues);
    setParameters!(block.block, originalValues);
    return jacobian;
end
