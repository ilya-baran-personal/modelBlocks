# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks

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

outputs = Variables([
    Variable("x", 0.0, (0, 100), "?", ""),
    Variable("y", 0.0, (0, 100), "?", ""),
]);

blockWithOutputs = BlockWithOutputs(block, outputs, (variables, parameters, timeRange, solution, outputs) -> begin
        outputs.x = solution.X[end];
        outputs.y = solution.Y[end];
        return outputs;
    end);

@time vpop = generateVPop(blockWithOutputs, 0:20, [
    ("p1", 0., 5.),
    ("p2", 0., 5.)
], Dict(
    "x" => (1., 0.02),
    "y" => (0.2, 0.02)
), 5; MaxTime = 10);
