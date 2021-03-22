# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using OrdinaryDiffEq

struct Block
    variables::Variables
    parameters::Variables
    reactions::Array{AbstractReaction}
end

function computeDerivatives(block::Block, t::Number, x::Vector)::Vector
    derivatives = deepcopy(block.variables);
    variables = deepcopy(block.variables);
    copyto!(derivatives.values, zero(derivatives.values));
    copyto!(variables.values, x);
    for reaction = block.reactions
        apply!(reaction, derivatives, variables, block.parameters)
    end
    #println("Time = $t, x = $x, result = $(derivatives.values)");
    return derivatives.values;
end

function runBlock(block::Block, timeInterval::Tuple{Number, Number})
    floatInterval::Tuple{Float64, Float64} = convert(Tuple{Float64, Float64}, timeInterval);
    problem = ODEProblem((x, p, t) -> computeDerivatives(block, t, x), block.variables.values, floatInterval);
    solve(problem, AutoTsit5(Rosenbrock23()));
end
