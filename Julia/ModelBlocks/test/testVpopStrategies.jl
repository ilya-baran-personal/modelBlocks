# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks
using Test
using LinearAlgebra
using Plots

variables = Variables([
    Variable("X", 1, (-100, 100), "", ""),
    Variable("Y", 1, (-100, 100), "", ""),
]);

parameters = Variables([
    Variable("x", 0.25, (-2, 2), "", ""),
    Variable("y", 0.375, (-2, 2), "", ""),
]);

reactions = [
    GeneralReaction("Sines", (dvdt, t, v, p, e) -> begin
        dvdt.X += (p.x - v.X);
        dvdt.Y += (p.y - v.Y);
    end),
];

block = Block(variables, parameters, reactions);
setTimeRange!(block, 0:10);

outputs = Variables([
    Variable("o1", 0.0, (0, 100), "?", ""),
    Variable("o2", 0.0, (0, 100), "?", ""),
]);

setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
        outputs.o1 = solution.X[end] ^ 2 + solution.Y[end] ^ 2;
        outputs.o2 = 0;
        return outputs;
    end);

@time ppop = generatePPopFarthest([block], [
    ("x", -2., 2.),
    ("y", -2., 2.)
], Dict(
    "o1" => (0., 1.),
    "o2" => (0., 1.)
), 100; MaxTime = 120, DistanceFactor = 1.);

Plots.plot(ppop[1,:], ppop[2,:], seriestype = :scatter)

squareSums = ones(size(ppop, 2), 1) * sum(ppop .^ 2, dims=1);
distances = squareSums + squareSums' - 2 * ppop' * ppop;
distances = reshape(distances, 1, :)[:];

Plots.histogram(log.(distances[distances .> 1e-5]), xlims=(-8, 2), ylims=(0,1000))
