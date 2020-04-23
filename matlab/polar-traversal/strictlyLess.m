function [check] = strictlyLess(a, b, absEpsilon, relEpsilon)
    check = a < b && ~approximatelyEqual(a,b,absEpsilon,relEpsilon);
end