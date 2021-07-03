# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks
using Test
using LinearAlgebra

@test nBallVolume([2, 3], 2) ≈ [4 * pi, 9 * pi];
@test nBallVolume(2, 3) ≈ 32 * pi / 3;

variables = Variables([
    Variable("X", 1, (0, 100), "", ""),
    Variable("Y", 1, (0, 100), "", ""),
]);

parameters = Variables([
    Variable("p1", 3.31583, (0, 100), "", ""),
    Variable("p2", 0.912226, (0, 100), "", ""),
]);

reactions = [
    GeneralReaction("Sines", (dvdt, t, v, p, e) -> begin
        dvdt.X += p.p1 * v.Y;
        dvdt.Y -= p.p2 * v.X;
    end),
];

block = Block(variables, parameters, reactions);
setTimeRange!(block, 0:20);

outputs = Variables([
    Variable("x", 0.0, (0, 100), "?", ""),
    Variable("y", 0.0, (0, 100), "?", ""),
]);

setOutputDefinition!(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
        outputs.x = solution.X[end];
        outputs.y = solution.Y[end];
        return outputs;
    end);

@time ppop = generatePPop(block, [
    ("p1", 0., 5.),
    ("p2", 0., 5.)
], Dict(
    "x" => (1., 0.12),
    "y" => (0.2, 0.12)
), 8; MaxTime = 50);

ppop = expandPPop([block], [
    ("p1", 0., 5.),
    ("p2", 0., 5.)
], Dict(
    "x" => (1., 0.12),
    "y" => (0.2, 0.12)
), ppop, 0.3, 1000);

@time samples = subsamplePPop(ppop, block, ["p1", "p2"], fill(0.5, 2), [[1,1] [0,1]] / 100, 2000);

ppop2 = rand(2, 20000);
@time samples = subsamplePPop(ppop2, fill(0.5, 2), [[1,1] [0,1]] / 100, 2000);
