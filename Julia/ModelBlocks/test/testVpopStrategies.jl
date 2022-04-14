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
        outputs.o1 = cos(4 * sqrt(solution.X[end] ^ 2 + solution.Y[end] ^ 2));
        outputs.o2 = 0;
        return outputs;
    end);

@time ppop = generatePPopFarthest([block], [
    ("x", -2., 2.),
    ("y", -2., 2.)
], Dict(
    "o1" => (0.5, 0.5),
    "o2" => (0.5, 0.5)
), 100; MaxTime = 120, DistanceFactor = 0.5);

plt = plot(ppop[1,:], ppop[2,:], seriestype = :scatter, legend = false, aspect_ratio = :equal, size = (600, 600));
plot!(pi/8*cos.(0:0.1:(2*pi+0.1)), pi/8*sin.(0:0.1:(2*pi+0.1)));
plot!(3*pi/8*cos.(0:0.1:(2*pi+0.1)), 3*pi/8*sin.(0:0.1:(2*pi+0.1)));
plot!(5*pi/8*cos.(0:0.1:(2*pi+0.1)), 5*pi/8*sin.(0:0.1:(2*pi+0.1)));
plot!(max.(-2., min.(2., 7*pi/8*cos.(0:0.01:(2*pi+0.1)))), max.(-2., min.(2., 7*pi/8*sin.(0:0.01:(2*pi+0.1)))));
display(plt);

# squareSums = ones(size(ppop, 2), 1) * sum(ppop .^ 2, dims=1);
# distances = squareSums + squareSums' - 2 * ppop' * ppop;
# distances = reshape(distances, 1, :)[:];

# histogram(log.(distances[distances .> 1e-5]), xlims=(-8, 2), ylims=(0,1000))


# plot volume-growth curve
radii = 0:0.01:0.5;
areas = zeros(length(radii));
samples = 100000; # if checking outputs, you'll need fewer samples

n = size(ppop)[2];

for i = 1:samples
    s = rand(2) * 4 .- 2; # The 4 and -2 should become range widths and minima
# Check outputs here
    d = sum((ppop - repeat(s, 1, n)) .^ 2, dims = 1);
    closest = sqrt(minimum(d));
    areas[radii .> closest] .+= 4 / samples;
end

plt = plot(radii, areas, legend = false);
display(plt);