function [check] = approximatelyEqual(a, b, absEpsilon, relEpsilon)
% https://www.learncpp.com/cpp-tutorial/relational-operators-and-floating-point-comparisons/
    diff = abs(a-b);
    if diff <= absEpsilon
        check = true;
        return
    end
    check = (diff <= max(abs(a),abs(b))*relEpsilon);
end