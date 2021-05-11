# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.
module ModelBlocks

export Variable, Variables, variablesUnion, variablesSubtract, variablesToMatlab
include("variables.jl")

export AbstractReaction, SimpleReaction, GeneralRateReaction, GeneralReaction
include("reactions.jl")

export AbstractBlock, Block, getVariables, getParameters, runBlock, BlockWithBindings, BlockWithOutputs, getOutputs
export solutionToMatlab, solutionToMatrix, solutionToVariables, setParameter!, setParameters!
include("block.jl")

export localSensitivity
include("sensitivity.jl")

export fitParameters
include("fitting.jl")

export generateVPop
include("vpop.jl")

end
