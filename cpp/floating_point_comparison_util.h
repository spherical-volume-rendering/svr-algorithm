#ifndef SVR_FLOATING_POINT_COMPARISON_UTIL_H
#define SVR_FLOATING_POINT_COMPARISON_UTIL_H

#include "vec3.h"
#include <cmath>
#include <algorithm>

// A number of floating point comparison algorithms written for the spherical volume rendering algorithm.

namespace svr {
    // Epsilons used for floating point comparisons in Knuth's algorithm.
    constexpr double ABS_EPSILON = 1e-12;
    constexpr double REL_EPSILON = 1e-8;

    // Determines equality between two floating point numbers using a defaulted absolute and relative epsilon.
    // Related Boost document:
    //        https://www.boost.org/doc/libs/1_61_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/
    //        floating_point/floating_points_comparison_theory.html#equ1
    // Related reading:
    //        Donald. E. Knuth, 1998, Addison-Wesley Longman, Inc., ISBN 0-201-89684-2, Addison-Wesley Professional;
    //        3rd edition. (The relevant equations are in §4.2.2, Eq. 36 and 37.)
    inline bool isEqual(double a, double b) noexcept {
        const double diff = std::abs(a - b);
        if (diff <= ABS_EPSILON) { return true; }
        return diff <= std::max(std::abs(a), std::abs(b)) * REL_EPSILON;
    }

    // Overloaded version that checks for Knuth equality with vector cartesian coordinates.
    inline bool isEqual(const Vec3 &a, const Vec3 &b) noexcept {
        const double diff_x = std::abs(a.x() - b.x());
        const double diff_y = std::abs(a.y() - b.y());
        const double diff_z = std::abs(a.z() - b.z());
        if (diff_x <= ABS_EPSILON && diff_y <= ABS_EPSILON && diff_z <= ABS_EPSILON) { return true; }
        return diff_x <= std::max(std::abs(a.x()), std::abs(b.x())) * REL_EPSILON &&
               diff_y <= std::max(std::abs(a.y()), std::abs(b.y())) * REL_EPSILON &&
               diff_z <= std::max(std::abs(a.z()), std::abs(b.z())) * REL_EPSILON;
    }

    // Checks to see if a is strictly less than b using Knuth's algorithm.
    inline bool lessThan(double a, double b) noexcept {
        return a < b && !isEqual(a, b);
    }

} // svr

#endif //SVR_FLOATING_POINT_COMPARISON_UTIL_H
