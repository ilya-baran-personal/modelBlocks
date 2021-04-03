# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using OrdinaryDiffEq

struct Block
    variables::Variables
    parameters::Variables
    reactions::Array{AbstractReaction}
    Block(v,p,r) = new(deepcopy(v), deepcopy(p), deepcopy(r));
end

function computeDerivatives(block::Block, t::Number, x::Vector)::Vector
    derivatives = Variables(block.variables, zero(x));
    variables = Variables(block.variables, x);
    for reaction = block.reactions
        apply!(reaction, derivatives, t, variables, block.parameters);
    end
    #println("Time = $t, x = $x, result = $(derivatives.values)");
    return derivatives.values;
end

function runBlock(block::Block, timeRange::AbstractRange)
    floatInterval::Tuple{Float64, Float64} = convert(Tuple{Float64, Float64}, (minimum(timeRange), maximum(timeRange)));
    problem = ODEProblem((x, p, t) -> computeDerivatives(block, t, x), block.variables.values, floatInterval);
    solution = solve(problem, AutoTsit5(Rosenbrock23()));
end

function solutionToVariables(block::Block, timeRange::AbstractRange, solution)
    resultArray = Array{Array{Float64, 1}, 1}(undef, length(block.variables.variables));
    for index = 1:length(resultArray)
        resultArray[index] = solution(timeRange, idxs=index).u;
    end
    return Variables(block.variables, resultArray);
end

struct BlockWithOutputs
    block::Block
    outputs::Variables
    computeOutputs::Function # Takes variables, solution (in Variables form), and outputs, and returns outputs
end

function getOutputs(blockWithOutputs::BlockWithOutputs, timeRange::AbstractRange)
    solution = runBlock(blockWithOutputs.block, timeRange);
    return blockWithOutputs.computeOutputs(blockWithOutputs.block.variables, timeRange,
                                           solutionToVariables(blockWithOutputs.block, timeRange, solution),
                                           blockWithOutputs.outputs);
end