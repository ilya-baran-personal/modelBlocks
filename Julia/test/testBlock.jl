# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

include("../src/ModelBlocks.jl");

using .ModelBlocks
using Test
using Plots
using OrdinaryDiffEq

println("Block Test starting");

variables = Variables([
    Variable("X", 1, (0, 100), "kg", ""),
    Variable("Y", 1, (0, 100), "kg", ""),
    Variable("Z", 1, (0, 100), "kg", ""),
]);

parameters = Variables([
    Variable("Xin", 1.0, (0, 100), "kg/s", ""),
    Variable("Yin", 1.0, (0, 100), "kg/s", ""),
    Variable("Zout", 0.5, (0, 100), "kg/s", ""),
    Variable("Kxy", 1.0, (0, 100), "?", ""),
]);

reactions = [
    SimpleReaction("X uptake", [], ["X"], "Xin"),
    SimpleReaction("Y uptake", [], ["Y"], "Yin"),
    SimpleReaction("Z clearance", ["Z"], [], "Zout"),
    SimpleReaction("Transformation", ["X", "Y"], ["Z"], "Kxy"),
];

block = Block(variables, parameters, reactions);

@time solution = runBlock(block, (0, 100));

println("Block Test finished");
