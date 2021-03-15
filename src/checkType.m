... Copyright (C) 2021 Ilya Baran.  This program is distributed under the
... terms of the MIT license.
function checkType(value, type, name)
if (~isa(value, type))
    error ('Expected %s to be a %s, got %s', name, type, class(value));
end
end

