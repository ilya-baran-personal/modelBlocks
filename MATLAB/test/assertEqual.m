... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
function assertEqual(expected, actual)
%ASSERTEQUAL asserts that two inputs are the same
if ~isequal(expected, actual)
    expectedString = evalc('disp(expected)');
    actualString = evalc('disp(actual)');
    error('Expected:\n%s\nactual:\n%s', expectedString, actualString);
end
end
