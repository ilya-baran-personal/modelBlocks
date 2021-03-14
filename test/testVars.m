addpath('../src');

varList = {};
varList = var('cats', 3, 0, 10, 'unitless', 'Number of cats').addTo(varList);
varList = var('dogs', 2, 0, 10, 'unitless', 'Number of dogs').addTo(varList);
varList = var('catFood', 3, 0, 10, 'kilograms').addTo(varList);

myVars = vars(varList)

assert(myVars.cats == 3);
assert(myVars.dogs == 2);
assert(myVars.values(2) == 2);
assert(myVars.variables{2}.value == 2);
assert(strcmp(myVars.catFood.units, 'kilograms'));
myVars.catFood = 2.5;
assert(myVars.catFood == 2.5);
