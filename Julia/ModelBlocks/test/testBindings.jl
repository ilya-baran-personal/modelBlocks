# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks
using Test
using OrdinaryDiffEq
using Plots

println("Bindings Test starting");

variables = Variables([
    Variable("X", 1, (0, 100), "", ""),
    Variable("Y", 1, (0, 100), "", ""),
]);

parameters = Variables([
    Variable("p1", 1.0, (0, 100), "", ""),
    Variable("p2", 2.0, (0, 100), "", ""),
    Variable("p3", 3.0, (0, 100), "", ""),
    Variable("p4", 4.0, (0, 100), "", ""),
]);

reactions = [
    GeneralReaction("Really random stuff", (dvdt, t, v, p, e) -> begin
        dvdt.X += p.p1 * v.Y;
        dvdt.Y -= p.p2 * v.X;
        @test p.p1 == 3;
        @test p.p3 == t / 5.;
    end),
];

block = Block(variables, parameters, reactions);

boundBlock = BlockWithBindings(block, Dict([
    ("p1", 3),
    ("p3", (t, p, v) -> t / 5.)
]));

@test getVariables(boundBlock) == getVariables(block);
@test getParameters(boundBlock).values == [2.0, 4.0];

@time solution = runBlock(boundBlock, 0:100);

println("Bindings Test finished");
