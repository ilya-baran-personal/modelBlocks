# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.
module ModelBlocks

export Variable, Variables, variablesUnion, variablesSubtract, renameVariables, variablesToMatlab
include("variables.jl")

export AbstractReaction, SimpleReaction, GeneralRateReaction, GeneralReaction
include("reactions.jl")

export AbstractBlock, Block, getVariables, getParameters, runBlock, BlockWithBindings, BlockOutputDefinition, BlockExtraData, computeOutputs
export solutionToMatlab, solutionToMatrix, solutionToVariables, setParameter!, setParameters!, getExtraData, getTimeRange, setTimeRange!
export getOutputDefinition, setOutputDefinition!, getDiscontinuities, setDiscontinuities!, BlockCombo, renameOutputs!
include("block.jl")

export localSensitivity
include("sensitivity.jl")

export fitParameters, fitParameters!
include("fitting.jl")

export generatePPop, subsamplePPop, nBallVolume, normalPDF, weightedSample
include("vpop.jl")

end
