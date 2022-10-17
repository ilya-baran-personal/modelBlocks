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
    # Variable("z", 0.375, (-2, 2), "", ""),
    # Variable("w", 0.375, (-2, 2), "", ""),
    # Variable("q", 0.375, (-2, 2), "", ""),
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

r = 3.91;

setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
        outputs.o1 = cos(r * sqrt(solution.X[end] ^ 2 + solution.Y[end] ^ 2));
        outputs.o2 = 0.5;
        return outputs;
    end);

parameterBounds =  [
    ("x", -2., 2.),
    ("y", -2., 2.),
    # ("z", 1., 2.),
    # ("w", 1., 2.),
    # ("q", 1., 2.),
];

outputBounds = Dict(
    "o1" => (0.5, 0.5),
    "o2" => (0.5, 0.5)
);

@time ppop = generatePPop(block, parameterBounds, outputBounds, 100; MaxTime = 60);
@time ppopFarthest = generatePPopFarthest([block], parameterBounds, outputBounds, 100; MaxTime = 60, DistanceFactor = 0.5, threads = 1);
@time ppopFarthestMT = generatePPopFarthest([block], parameterBounds, outputBounds, 100; MaxTime = 30, DistanceFactor = 0.5, threads = 4);
@time ppopFarthestMTNoPerturb = generatePPopFarthest([block], parameterBounds, outputBounds, 100; MaxTime = 30, DistanceFactor = 0.5, threads = 4, perturb = false);

for pts in (ppopFarthest, ppop)
    plt = plot(pts[1,:], pts[2,:], seriestype = :scatter, legend = false, aspect_ratio = :equal, size = (600, 600));
    plot!(pi/(2*r)*cos.(0:0.1:(2*pi+0.1)), pi/(2*r)*sin.(0:0.1:(2*pi+0.1)));
    plot!(3*pi/(2*r)*cos.(0:0.1:(2*pi+0.1)), 3*pi/(2*r)*sin.(0:0.1:(2*pi+0.1)));
    plot!(5*pi/(2*r)*cos.(0:0.1:(2*pi+0.1)), 5*pi/(2*r)*sin.(0:0.1:(2*pi+0.1)));
    plot!(max.(-2., min.(2., 7*pi/(2*r)*cos.(0:0.01:(2*pi+0.1)))), max.(-2., min.(2., 7*pi/(2*r)*sin.(0:0.01:(2*pi+0.1)))));
    plt = plot(plt, thickness_scaling=2, tickfontsize=10/2, labelfontsize=14/2, colorbar_tickfontsize=8/2, reuse=false);
    display(plt);
end

# Add corners: 
# for i in 1:2 for j in 1:2 for k in 1:2 for q in 1:2 for r in 1:2 ppopFarthest = hcat(ppopFarthest, [i;j;k;q;r]) end end end end end

# plot volume-growth curve
(radii, areas) = computeDistanceCurve([block], parameterBounds, outputBounds, ppop, 0.15; samples = 1000);
(radii, areasFarthest) = computeDistanceCurve([block], parameterBounds, outputBounds, ppopFarthest, 0.15; samples = 1000);
(radii, areasFarthestMT) = computeDistanceCurve([block], parameterBounds, outputBounds, ppopFarthestMT, 0.15; samples = 10000);
(radii, areasFarthestMTNoPerturb) = computeDistanceCurve([block], parameterBounds, outputBounds, ppopFarthestMTNoPerturb, 0.15; samples = 10000);

plt = plot(radii, hcat(areas, areasFarthest), label = ["SA" "FP"]);
plt = plot(plt, thickness_scaling=2, tickfontsize=10/2, labelfontsize=14/2, colorbar_tickfontsize=8/2, reuse=false);
display(plt);

# for i in 130:180
#     (radii, areas) = computeDistanceCurve([block], parameterBounds, outputBounds, ppop[:,1:i], 0.15; samples = 3000);
#     println("i = $i diff = $(sum(areasFarthest) - sum(areas))");
# end