# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.
module ModelBlocks

export Variable, Variables, variablesUnion
include("variables.jl")

export AbstractReaction, SimpleReaction, GeneralRateReaction, GeneralReaction
include("reactions.jl")

export Block, getVariables, getParameters, runBlock, BlockWithBindings, BlockWithOutputs, getOutputs
export solutionToMatlab, solutionToMatrix, solutionToVariables
include("block.jl")

export localSensitivity
include("sensitivity.jl")

export fitParameters
include("fitting.jl")

end
