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
    GeneralReaction("Really random stuff", (dvdt, t, v, p) -> begin
        dvdt.X += 0.1;
        dvdt.Y -= v.Z * 0.1;
    end)
];

block = Block(variables, parameters, reactions);

@time solution = runBlock(block, 0:100);

println("Block Test finished");
