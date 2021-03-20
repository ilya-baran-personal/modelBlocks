... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
clear variables;
addpath('../src');

% Define variables
varList = {};
varList = var('water', 100, 0, 1000, 'kilograms', 'Water in the pool').addTo(varList);
varList = var('salt', 0, 0, 10, 'kilograms', 'Salt in the pool').addTo(varList);

variables = vars(varList);

% Define parameters
varList = {};
varList = var('water_in', 0.5, 0, 10, 'kilograms/second', 'Water added').addTo(varList);
varList = var('salt_in', 0.1, 0, 10, 'kilograms/second', 'Salt added').addTo(varList);
varList = var('out_ratio', 0.01, 0, 10, '1/second', ...
    'What fraction leaves pool every second').addTo(varList);

parameters = vars(varList);

b = block(variables, parameters);
% Define reactions
b = b.addReaction(basicReaction('water_in',  {},        {'water'}, 'water_in'));
b = b.addReaction(basicReaction('salt_in',   {},        {'salt'},  'salt_in'));
b = b.addReaction(basicReaction('water_out', {'water'}, {},        'out_ratio'));
b = b.addReaction(basicReaction('salt_out',  {'salt'},  {},        'out_ratio'));

tic;
result = b.run([0, 1000]); % run for 10s
toc
tic;
result = ode15s(@(t,y) -0.01 * y, [0, 1000], [20; 100]);
toc

