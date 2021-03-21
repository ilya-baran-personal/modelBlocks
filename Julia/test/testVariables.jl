# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

include("../src/ModelBlocks.jl");

using .ModelBlocks
using Test

println("Variables Test starting");

variableArray = [
    Variable("foo", 3, (2, 4), "meters", " Basic foo variable"),
    Variable("bar", 4, (2, 5), "meters", " Basic bar variable"),
    Variable("m_3", 2, (1, 4), "moles", " This is m_3"),
];
variables = Variables(variableArray);

@test variables.foo == 3;
@test variables.bar == 4;
@test variables.m_3 == 2;
@test variables.values == [3, 4, 2];
@test variables.nameToIndex["foo"] == 1;

println("Variables Test finished");