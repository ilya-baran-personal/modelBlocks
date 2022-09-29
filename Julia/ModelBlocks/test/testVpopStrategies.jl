# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks
using Test
using LinearAlgebra
using Plots

#plotly();

variables = Variables([
    Variable("X", 1, (-100, 100), "", ""),
    Variable("Y", 1, (-100, 100), "", ""),
]);

parameters = Variables([
    Variable("x", 0.25, (-2, 2), "", ""),
    Variable("y", 0.375, (-2, 2), "", ""),
    Variable("z", 0.375, (-2, 2), "", ""),
    Variable("w", 0.375, (-2, 2), "", ""),
    Variable("q", 0.375, (-2, 2), "", ""),
]);

reactions = [
    # GeneralReaction("Sines", (dvdt, t, v, p, e) -> begin
    #     dvdt.X += (p.x - v.X);
    #     dvdt.Y += (p.y - v.Y);
    # end),
];

block = Block(variables, parameters, reactions);
setTimeRange!(block, 0:0);

outputs = Variables([
    Variable("o1", 0.0, (0, 100), "?", ""),
    Variable("o2", 0.0, (0, 100), "?", ""),
]);

setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
        outputs.o1 = 0.5;#cos(4 * sqrt(solution.X[end] ^ 2 + solution.Y[end] ^ 2));
        outputs.o2 = 0.5;
        return outputs;
    end);

parameterBounds =  [
    ("x", 1., 2.),
    ("y", 1., 2.),
    ("z", 1., 2.),
    ("w", 1., 2.),
    ("q", 1., 2.),
];

outputBounds = Dict(
    "o1" => (0.5, 0.5),
    "o2" => (0.5, 0.5)
);

@time ppopFarthest = generatePPopFarthest([block], parameterBounds, outputBounds, 1000; MaxTime = 120, DistanceFactor = 0.5, threads = 1);
@time ppop = generatePPop(block, parameterBounds, outputBounds, 1000; MaxTime = 120);

plt = plot(ppopFarthest[1,:], ppopFarthest[2,:], seriestype = :scatter, legend = false, aspect_ratio = :equal, size = (600, 600));
# plot!(pi/8*cos.(0:0.1:(2*pi+0.1)), pi/8*sin.(0:0.1:(2*pi+0.1)));
# plot!(3*pi/8*cos.(0:0.1:(2*pi+0.1)), 3*pi/8*sin.(0:0.1:(2*pi+0.1)));
# plot!(5*pi/8*cos.(0:0.1:(2*pi+0.1)), 5*pi/8*sin.(0:0.1:(2*pi+0.1)));
# plot!(max.(-2., min.(2., 7*pi/8*cos.(0:0.01:(2*pi+0.1)))), max.(-2., min.(2., 7*pi/8*sin.(0:0.01:(2*pi+0.1)))));
display(plt);

# Add corners: 
# for i in 1:2 for j in 1:2 for k in 1:2 for q in 1:2 for r in 1:2 ppopFarthest = hcat(ppopFarthest, [i;j;k;q;r]) end end end end end

# plot volume-growth curve
(radii, areas) = computeDistanceCurve([block], parameterBounds, outputBounds, ppop, 0.3; samples = 30000);
(radii, areasFarthest) = computeDistanceCurve([block], parameterBounds, outputBounds, ppopFarthest, 0.3; samples = 30000);

plt = plot(radii, hcat(areas, areasFarthest));
display(plt);
