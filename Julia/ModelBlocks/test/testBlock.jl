# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks
using Test
using OrdinaryDiffEq
using Plots

println("Block Test starting");

variables = Variables([
    Variable("X", 1, (0, 100), "kg", ""),
    Variable("Y", 1, (0, 100), "kg", ""),
    Variable("Z", 1, (0, 100), "kg", ""),
]);

parameters = Variables([
    Variable("Xin", 1.5, (0, 100), "kg/s", ""),
    Variable("Yin", 1.0, (0, 100), "kg/s", ""),
    Variable("Zout", 1.0, (0, 100), "kg/s", ""),
    Variable("Xout", .2, (0, 100), "kg/s", ""),
    Variable("Kxy", 0.5, (0, 100), "?", ""),
    Variable("R", .1, (0, 100), "?", ""),
]);

reactions = [
    SimpleReaction("X uptake", [], ["X"], "Xin"),
    SimpleReaction("Y uptake", [], ["Y"], "Yin"),
    SimpleReaction("Z clearance", ["Z"], [], "Zout"),
    SimpleReaction("X clearance", ["X"], [], "Xout"),
    SimpleReaction("Transformation", ["X", "Y"], ["Z"], "Kxy"),
    GeneralRateReaction("Random stuff", ["X"], ["Y"], (t, v, p) -> v.X * cos(t) * p.R),
    GeneralReaction("Really random stuff", (dvdt, t, v, p, e) -> begin
        dvdt.X += 0.1;
        dvdt.Y -= v.Z * 0.1;
        e["Extra"] = 3;
    end),
    GeneralRateReaction("Random stuff", ["X"], ["Y"], (t, v, p, e) -> e["Extra"]),
];

block = Block(variables, parameters, reactions);
setTimeRange!(block, 0:100);

@time solution = runBlock(block);

outputs = Variables([
    Variable("vector", 0.0, (0, 100), "?", ""),
    Variable("scalar", 0.0, (0, 100), "?", ""),
]);

setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
    outputs.vector = [solution.X[1], solution.X[5]];
    outputs.scalar = solution.Y[end];

    return outputs;
end);

@time fit = fitParameters([block], [
    ("Xin", 0., 10.0),
    ("Yin", 0., 10.0),
], Dict(
    "vector" => ([1, 2], [1, 1]),
    "scalar" => (3, 1)
); MaxTime = 3);

println("Block Test finished");
