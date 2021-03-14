... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
clear variables;
addpath('../src');

varList = {};
varList = var('cats', 3, 0, 10, 'unitless', 'Number of cats').addTo(varList);
varList = var('dogs', 2, 0, 10, 'unitless', 'Number of dogs').addTo(varList);
varList = var('catFood', 3, 0, 10, 'kilograms').addTo(varList);

myVars = vars(varList)

assertEqual(3, myVars.cats);
assertEqual(2, myVars.dogs);
assertEqual([3, 2, 3], myVars.values);
assertEqual(2, myVars.variables{2}.value);
assertEqual('kilograms', myVars.catFood.units);
myVars.catFood = 2.5;
assertEqual(2.5, myVars.catFood);

% Test invalid name
success = true;
try
    var('3cats', 1, 0, 1, 'units') % Invalid name
    success = false;
catch e % expected
end
assert(success, 'var should have complained about invalid name');

% Test duplicate variables in vars
varList = var('cats', 3, 0, 10, 'unitless', 'Number of cats').addTo(varList);
try
    vars(varList)
    success = false;
catch e % expected
end
assert(success, 'vars should have complained about duplicate variables');
