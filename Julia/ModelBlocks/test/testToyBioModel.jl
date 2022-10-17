# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks
using Test
using LinearAlgebra
using Interpolations
using Plots


# Part 1: define model - variables, parameters, reactions, block, setTimeRange, outputs

println("Toy model test starting");

variables = Variables([
    Variable("X", 4., (0, 100), "mol/L", ""),
    Variable("Y", 0.0, (0, 100), "mol/L", "Enzyme that metabolizes X"),
]);

# parameters = Variables([
#     Variable("k_synX", 0.1, (0, 1000), "mol/L/s", "Synthesis rate of X"),
#     Variable("k_metX", 0.0001, (0, 1000), "1/s", "Metabolism rate of X"),
#     Variable("k_synY", 0.01, (0, 1000), "mol/L/s", "Synthesis rate of Y"),
#     Variable("k_metY", .01, (0, 1000), "1/s", "Metabolism rate of Y"),
#     Variable("k_elXY", .005, (0, 1000), "mol/L/s", "Elimination rate of X by Y"),
#     Variable("k_stimY", 0.09, (0, 1000), "mol/L/s", "The rate at which X stimulates production of Y"),
# ]);


parameters = Variables([
    Variable("k_synX", 0.1, (0, 1000), "mol/L/s", "Synthesis rate of X"),
    Variable("k_metX", 0.1, (0, 1000), "1/s", "Metabolism rate of X"),
    Variable("k_synY", 0.1, (0, 1000), "mol/L/s", "Synthesis rate of Y"),
    Variable("k_metY", .05, (0, 1000), "1/s", "Metabolism rate of Y"),
    Variable("k_elXY", .1, (0, 1000), "mol/L/s", "Elimination rate of X by Y"),
    Variable("k_stimY", 0.9, (0, 1000), "mol/L/s", "The rate at which X stimulates production of Y"),
    Variable("k_xyy", 0.0, (0, 1000), "mol/L/s", "The rate at which X stimulates production of Y"),
]);

reactions = [
    GeneralReaction("ToyModel", (dvdt, t, v, p, e) -> begin
        dvdt.X = p.k_synX - p.k_metX * v.X - p.k_elXY * v.X * v.Y + p.k_xyy * v.Y^2;
        dvdt.Y = p.k_synY - p.k_metY * v.Y + p.k_stimY * v.X;

        dvdt.X += 10 * max(0., sin(t)) ^ 10.; # Every 2pi seconds, add an infusion of X (e.g. via food intake)
    end),
];

block = Block(variables, parameters, reactions);
setTimeRange!(block, 0:.05:500);

outputs = Variables([
    Variable("minX", 0.0, (0, 100), "?", ""),
    Variable("maxX", 0.0, (0, 100), "?", ""),
    Variable("maxY", 0.0, (0, 100), "?", ""),
]);

setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
    outputs.minX = minimum(solution.X);
    outputs.maxX = maximum(solution.X);
    outputs.maxY = maximum(solution.Y);

    return outputs;
end);

println("Toy model Test defined");

# compute outputs
@time outputs = computeOutputs(block)

solution = runBlock(block);
plot(solution)

# Part 2: generate PPop

parameterBounds = [
    # Ranges are min to max
    ("k_synX", 0.02, 0.5),
    ("k_metX", 0.02, 0.5),
    ("k_synY", 0.02, 0.5),
    ("k_metY", 0.02, 0.5),
    ("k_elXY", 0.01, 2.),
    ("k_stimY", 0.01, 2.),
    ("k_xyy", 0., 0.03),
];

outputBounds = Dict(
    # All outputs ranges for PPs presented as (MEAN, STD)
     "minX" => (3., 2.9),
     "maxX" => (15., 10.),
     "maxY" => (15., 10.),
 );

# Plotting interactive plots with PyPlot:
# pygui(true)
# maybe: ion()
# PyPlot.scatter(ppop[1,:],ppop[2,:])
# scatter3D(ppop[1,:],ppop[2,:], ppop[3,:])

println("Generate PPops using Farthest Point optimization");
@time ppopFarthest = generatePPopFarthest([block], parameterBounds, outputBounds, 100;
                                          MaxTime = 120, DistanceFactor = 0.1, threads = 8); # num of pts, total max time for all PPs generation, distance to nearest existing pt is weighted relative to staying wihtin the output bounds

plt = plot(ppopFarthest[1,:], ppopFarthest[2,:], seriestype = :scatter, legend = false, size = (600, 600));

# # log plot
# #plt = plot(log.(ppopFarthest[1,:]), log.(ppopFarthest[2,:]), seriestype = :scatter, legend = false, aspect_ratio = :equal, size = (600, 600));
display(plt);

println("Generate PPops");

@time ppop = generatePPop(block, parameterBounds, outputBounds, 100; MaxTime = 120);

plt = plot(ppop[1,:], ppop[2,:], seriestype = :scatter, legend = false, size = (600, 600));

# # log plot
# #plt = plot(log.(ppop[1,:]), log.(ppop[2,:]), seriestype = :scatter, legend = false, aspect_ratio = :equal, size = (600, 600));
# display(plt);

# (radii, areas) = computeDistanceCurve([block], parameterBounds, outputBounds, ppop, .1; samples = 10000);
# (radii, areasFarthest) = computeDistanceCurve([block], parameterBounds, outputBounds, ppopFarthest, .1; samples = 10000);

# plt = plot(radii, hcat(areas, areasFarthest), legend = false, label = ["SA" "FP"]);
# plt = plot(plt, thickness_scaling=2, tickfontsize=10/2, labelfontsize=14/2, colorbar_tickfontsize=8/2, reuse=false);
# display(plt);

# # Print the integral of the difference
# display((sum(areasFarthest) - sum(areas)) * (radii[2] - radii[1]));


# for i in 120:180
#     (radii, areas) = computeDistanceCurve([block], parameterBounds, outputBounds, ppop[:,1:i], 0.15; samples = 3000);
#     println("i = $i diff = $(sum(areasFarthest) - sum(areas))");
# end