# Copyright (C) 2021 Ilya Baran.  This program is distributed under the terms of the MIT license.

using ModelBlocks
using Test

println("Variables Test starting");

variableArray = [
    Variable("foo", 3, (2, 4), "meters", "Basic foo variable"),
    Variable("bar", 4, (2, 5), "meters", "Basic bar variable"),
    Variable("m_3", 2, (1, 4), "moles", "This is m_3"),
];
variables = Variables(variableArray);

@test variables.foo == 3;
@test variables.bar == 4;
@test variables.m_3 == 2;
@test variables.values == [3, 4, 2];
@test variables.nameToIndex["foo"] == 1;
@test variables["foo"] == 3;

copyto!(variables.values, [4, 2, 1]);
@test variables.values == [4, 2, 1];

@test_throws ErrorException Variable("3foo", 3, (2, 4), "meters", "Invalid identifier")
@test_throws ErrorException Variable("foo", 3, (5, 4), "meters", "Invalid range")

resize!(variableArray, 4);
variableArray[4] = Variable("foo", 3, (2, 4), "meters", "Duplicate foo variable");
@test_throws ErrorException Variables(variableArray);

# Test combineVariables

v1 = Variables([
    Variable("a", 1, (1, 4), "", "A"),
    Variable("b", 2, (1, 5), "", "B"),
    Variable("c", 3, (1, 6), "", "C"),
]);

v2 = Variables([
    Variable("b", 5, (1, 16), "", "B"),
    Variable("d", 7, (1, 14), "", "D"),
    Variable("c", 6, (1, 15), "", "C"),
]);

vUnion = variablesUnion(v1, v2);
@test vUnion.values == [1, 5, 6, 7];

println("Variables Test finished");