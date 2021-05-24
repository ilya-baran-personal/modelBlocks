# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks
using Test

println("Combo Test starting");

variables1 = Variables([
    Variable("X", 1, (0, 100), "kg", ""),
]);

parameters1 = Variables([
    Variable("p1", 0.0, (0, 100), "kg/s", ""),
    Variable("p2", 0.0, (0, 100), "kg/s", ""),
    Variable("p3", 0.0, (0, 100), "kg/s", ""),
]);

reactions1 = [
    GeneralRateReaction("R", [], ["X"], (t, v, p) -> p.p1 + p.p2 * v.X + p.p3 * t),
];

block1 = Block(variables1, parameters1, reactions1);

variables2 = Variables([
    Variable("Y", 1, (0, 100), "kg", ""),
]);

parameters2 = Variables([
    Variable("p3", 0.0, (0, 100), "kg/s", ""),
    Variable("p4", 2.0, (0, 100), "kg/s", ""),
]);

reactions2 = [
    GeneralRateReaction("R", [], ["Y"], (t, v, p) -> mod(t, 10) < 0.05 ? p.p4 : 0),
];

block2 = Block(variables2, parameters2, reactions2);

extraParameters = Variables([
    Variable("p5", 1.0, (0, 100), "kg/s", ""),
]);

combo = BlockCombo(
    [("first", block1), ("second", block2)],
    [
        ("first", "p1", (t, v, p) -> v.Y - v.X * p.p5)
    ],
    extraParameters
);

setTimeRange!(combo, 0:100);
setDiscontinuities!(combo, (p, timeRange) -> vcat([[t * 10, t * 10 + 0.05] for t in 0:10]...));

@time solution = runBlock(combo);

println("Combo Test finished");