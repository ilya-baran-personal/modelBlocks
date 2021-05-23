# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

import FiniteDiff

function localSensitivity(block::AbstractBlock)
    originalValues = Vector{Float64}(getParameters(block).values);
    lambda = (x) -> begin
        setParameters!(block, x);
        outputs = computeOutputs(block);
        return Vector{Float64}(outputs.values);
    end
    jacobian = FiniteDiff.finite_difference_jacobian(lambda, originalValues);
    setParameters!(block, originalValues);
    return jacobian;
end
