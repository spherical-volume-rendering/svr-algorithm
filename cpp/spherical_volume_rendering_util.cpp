#include "spherical_volume_rendering_util.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

// Epsilons used for floating point comparisons in Knuth's algorithm.
const double ABS_EPSILON = 1e-12;
const double REL_EPSILON = 1e-8;

// The type corresponding to the voxel(s) with the minimum tMax value for a given traversal.
enum VoxelIntersectionType {
    None = 0,
    Radial = 1,
    Angular = 2,
    Azimuthal = 3,
    RadialAngular = 4,
    RadialAzimuthal = 5,
    AngularAzimuthal = 6,
    RadialAngularAzimuthal = 7
};

// The parameters returned by radialHit().
struct RadialHitParameters {
    // The time at which a hit occurs for the ray at the next point of intersection with a radial section.
    // This is always calculated from the ray origin.
    double tMaxR;

    // The voxel traversal value of a radial step: 0, +1, -1. This is added to the current radial voxel.
    int tStepR;

    // Determine whether the current voxel traversal was a 'transition'. This is necessary to determine when
    // The radial steps should go from negative to positive or vice-versa,
    // since the radial voxels go from 1..N, where N is the number of radial sections.
    bool previous_transition_flag;

    // Determines whether the current hit is within time bounds (t, t_end).
    bool within_bounds;

    // Determines whether the current voxel hit has caused a step outside of the spherical voxel grid.
    bool exits_voxel_bounds;
};

// The parameters returned by angularHit().
struct AngularHitParameters {
    // The time at which a hit occurs for the ray at the next point of intersection with an angular section.
    // This is always calculated from the ray origin.
    double tMaxTheta;

    // The voxel traversal value of an angular step: 0, +1, -1. This is added to the current angular voxel.
    int tStepTheta;

    // Determines whether the current hit is within time bounds (t, t_end).
    bool within_bounds;
};

// The parameters returned by azimuthalHit().
struct AzimuthalHitParameters {
    // The time at which a hit occurs for the ray at the next point of intersection with an azimuthal section.
    // This is always calculated from the ray origin.
    double tMaxPhi;

    // The voxel traversal value of an azimuthal step: 0, +1, -1. This is added to the current azimuthal voxel.
    int tStepPhi;

    // Determines whether the current hit is within time bounds (t, t_end).
    bool within_bounds;
};

// A generalized set of hit parameters.
struct GenHitParameters {
    double tMax;
    int tStep;
    bool within_bounds;
};

// Determines equality between two floating point numbers in two steps. First, it uses the absolute epsilon, then it
// uses a modified version of an algorithm developed by Donald Knuth (which in turn relies upon relative epsilon).
// Provides default values for the absolute and relative epsilon. The "Kn" in the function name is short for Knuth.
// Related Boost document:
//        https://www.boost.org/doc/libs/1_61_0/libs/test/doc/html/boost_test/testing_tools/extended_comparison/
//        floating_point/floating_points_comparison_theory.html#equ1
// Related reading:
//        Donald. E. Knuth, 1998, Addison-Wesley Longman, Inc., ISBN 0-201-89684-2, Addison-Wesley Professional;
//        3rd edition. (The relevant equations are in §4.2.2, Eq. 36 and 37.)
inline bool isKnEqual(double a, double b) noexcept {
    const double diff = std::abs(a - b);
    if (diff <= ABS_EPSILON) { return true; }
    return diff <= std::max(std::abs(a), std::abs(b)) * REL_EPSILON;
}

// Overloaded version that checks for Knuth equality with vector cartesian coordinates.
inline bool isKnEqual(const Vec3 &a, const Vec3 &b) noexcept {
    const double diff_x = std::abs(a.x() - b.x());
    const double diff_y = std::abs(a.y() - b.y());
    const double diff_z = std::abs(a.z() - b.z());
    if (diff_x <= ABS_EPSILON && diff_y <= ABS_EPSILON && diff_z <= ABS_EPSILON) { return true; }
    return diff_x <= std::max(std::abs(a.x()), std::abs(b.x())) * REL_EPSILON &&
           diff_y <= std::max(std::abs(a.y()), std::abs(b.y())) * REL_EPSILON &&
           diff_z <= std::max(std::abs(a.z()), std::abs(b.z())) * REL_EPSILON;
}

// Uses the Knuth algorithm in KnEqual() to ensure that a is strictly less than b.
inline bool KnLessThan(double a, double b) noexcept {
    return a < b && !isKnEqual(a, b);
}

// Determines whether a radial hit occurs for the given ray. A radial hit is considered an intersection with
// the ray and a radial section. This follows closely the mathematics presented in:
// http://cas.xav.free.fr/Graphics%20Gems%204%20-%20Paul%20S.%20Heckbert.pdf
// Input:
//    ray: The given ray to check for intersection.
//    grid: The grid that the ray is intersecting with.
//    current_voxel_ID_r: The current radial voxel ID.
//    ray_sphere_vector_dot: The dot product of the vector difference between the ray origin and the sphere origin.
//    t: The current time.
//    t_end: The maximum allowable time before the minimum of (sphere, grid) exit.
//    v: The dot product between the ray's unit direction and the ray_sphere_vector.
//    prev_transition_flag: Determines whether the previous radial traversal was a 'transition'. A transition
//                          is defined as the change in sign of tStepR. Another way this can be determined is
//                          sequential hits with equal radii.
//
// Returns: The corresponding radial hit parameters.
RadialHitParameters radialHit(const Ray &ray, const SphericalVoxelGrid &grid, int current_voxel_ID_r,
                              double ray_sphere_vector_dot, double t, double t_end, double v,
                              bool previous_transition_flag) noexcept {
    const double current_radius = grid.sphereMaxRadius() - grid.deltaRadius() * (current_voxel_ID_r - 1);
    double r_a = std::max(current_radius - grid.deltaRadius(), grid.deltaRadius());

    double r_b;
    if (!previous_transition_flag) {
        // To find the next radius, we need to check the previous_transition_flag:
        // In the case that the ray has sequential hits with equal radii, e.g.
        // the innermost radial disc, this ensures that the proper radii are being checked.
        r_b = std::min(current_radius + grid.deltaRadius(), grid.sphereMaxRadius());
    } else { r_b = std::min(current_radius, grid.sphereMaxRadius()); }
    // Find the intersection times for the ray and the previous and next radial discs.
    const double ray_sphere_dot_minus_v_squared = ray_sphere_vector_dot - v * v;
    double discriminant_a = r_a * r_a - ray_sphere_dot_minus_v_squared;
    if (discriminant_a < 0.0) {
        r_a += grid.deltaRadius();
        discriminant_a = r_a * r_a - ray_sphere_dot_minus_v_squared;
    }
    const double d_a = std::sqrt(discriminant_a);
    std::vector<double> intersection_times(4);
    intersection_times[0] = ray.timeOfIntersectionAt(v - d_a);
    intersection_times[1] = ray.timeOfIntersectionAt(v + d_a);

    const double discriminant_b = r_b * r_b - ray_sphere_dot_minus_v_squared;
    if (discriminant_b >= 0.0) {
        const double d_b = std::sqrt(discriminant_b);
        intersection_times[2] = ray.timeOfIntersectionAt(v - d_b);
        intersection_times[3] = ray.timeOfIntersectionAt(v + d_b);
    }
    std::vector<double> times_gt_t = intersection_times;
    times_gt_t.erase(std::remove_if(times_gt_t.begin(), times_gt_t.end(), [t](double i) {
        return i < t || isKnEqual(i, t);
    }), times_gt_t.end());
    RadialHitParameters radial_params;
    bool t_within_bounds = false;
    if (times_gt_t.size() >= 2 && isKnEqual(intersection_times[0], intersection_times[1])) {
        // Ray is tangent to the circle, i.e. two intersection times are equal.
        radial_params.tMaxR = times_gt_t[0];
        const BoundVec3 p = ray.pointAtParameter(radial_params.tMaxR);
        const double p_x_new = p.x() - grid.sphereCenter().x();
        const double p_y_new = p.y() - grid.sphereCenter().y();
        const double p_z_new = p.z() - grid.sphereCenter().z();
        const double r_new = std::sqrt(p_x_new * p_x_new + p_y_new * p_y_new + p_z_new * p_z_new);

        radial_params.tStepR = 0;
        t_within_bounds = KnLessThan(t, radial_params.tMaxR) && KnLessThan(radial_params.tMaxR, t_end);
        if (isKnEqual(current_radius, r_new)) {
            radial_params.previous_transition_flag = true;
        } else {
            radial_params.previous_transition_flag = false;
        }
    } else if (times_gt_t.empty()) {
        // No intersection.
        radial_params.tMaxR = std::numeric_limits<double>::infinity();
        radial_params.tStepR = 0;
        radial_params.previous_transition_flag = false;
    } else {
        // Radial intersection.
        radial_params.tMaxR = times_gt_t[0];
        const BoundVec3 p = ray.pointAtParameter(radial_params.tMaxR);
        const double p_x_new = p.x() - grid.sphereCenter().x();
        const double p_y_new = p.y() - grid.sphereCenter().y();
        const double p_z_new = p.z() - grid.sphereCenter().z();
        const double r_new = std::sqrt(p_x_new * p_x_new + p_y_new * p_y_new + p_z_new * p_z_new);

        const bool radial_difference_is_near_zero = isKnEqual(r_new, current_radius);
        if (radial_difference_is_near_zero) {
            radial_params.previous_transition_flag = true;
        } else {
            radial_params.previous_transition_flag = false;
        }

        if (KnLessThan(r_new, current_radius) && !radial_difference_is_near_zero) {
            radial_params.tStepR = 1;
        } else {
            radial_params.tStepR = -1;
        }
        t_within_bounds = KnLessThan(t, radial_params.tMaxR) && KnLessThan(radial_params.tMaxR, t_end);
    }
    radial_params.exits_voxel_bounds = (current_voxel_ID_r + radial_params.tStepR) == 0;
    radial_params.within_bounds = t_within_bounds && !radial_params.exits_voxel_bounds;
    return radial_params;
}

// A generalized version of the latter half of the angular and azimuthal hit parameters. Since the only difference
// is the 2-d plane for which they exist in, this portion can be generalized to a single function call. The calculations
// presented below follow closely the works of [Foley et al, 1996], [O'Rourke, 1998].
// Quick reference: http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
GenHitParameters generalizedPlaneHit(const Ray &ray, double perp_uv_min, double perp_uv_max, double perp_uw_min,
                                     double perp_uw_max, double perp_vw_min, double perp_vw_max, const BoundVec3 &p,
                                     const FreeVec3 &v, double t, double t_end, double ray_plane_dir,
                                     std::size_t num_voxels) noexcept {
    const bool is_parallel_min = isKnEqual(perp_uv_min, 0.0);
    const bool is_parallel_max = isKnEqual(perp_uv_max, 0.0);
    double a, b;
    bool is_intersect_min, is_intersect_max;
    double t_min = 0.0;
    if (!is_parallel_min) {
        const double inv_perp_uv_min = 1.0 / perp_uv_min;
        a = perp_vw_min * inv_perp_uv_min;
        b = perp_uw_min * inv_perp_uv_min;
        if ((KnLessThan(a, 0.0) || KnLessThan(1.0, a)) ||
            KnLessThan(b, 0.0) || KnLessThan(1.0, b)) {
            is_intersect_min = false;
        } else {
            is_intersect_min = true;
            const BoundVec3 p_min(p.x() + v.x() * b, p.y() + v.y() * b, p.z() + v.z() * b);
            t_min = ray.timeOfIntersectionAt(p_min);
        }
    } else {
        is_intersect_min = false;
    }
    double t_max = 0.0;
    if (!is_parallel_max) {
        const double inv_perp_uv_max = 1.0 / perp_uv_max;
        a = perp_vw_max * inv_perp_uv_max;
        b = perp_uw_max * inv_perp_uv_max;
        if ((KnLessThan(a, 0.0) || KnLessThan(1.0, a)) ||
            KnLessThan(b, 0.0) || KnLessThan(1.0, b)) {
            is_intersect_max = false;
        } else {
            is_intersect_max = true;
            const BoundVec3 p_max(p.x() + v.x() * b, p.y() + v.y() * b, p.z() + v.z() * b);
            t_max = ray.timeOfIntersectionAt(p_max);
        }
    } else {
        is_intersect_max = false;
    }
    GenHitParameters params;
    if (is_intersect_max && !is_intersect_min && KnLessThan(t_max, t_end) && KnLessThan(t, t_max)) {
        params.tStep = 1;
        params.tMax = t_max;
        params.within_bounds = true;
        return params;
    }
    if (is_intersect_min && !is_intersect_max && KnLessThan(t_min, t_end) && KnLessThan(t, t_min)) {
        params.tStep = -1;
        params.tMax = t_min;
        params.within_bounds = true;
        return params;
    }
    if (is_intersect_max && is_intersect_min) {
        const bool t_min_within_bounds = KnLessThan(t, t_min) && KnLessThan(t_min, t_end);
        if (t_min_within_bounds && isKnEqual(t_min, t_max)) {

            params.tMax = t_max;
            if (!isKnEqual(ray_plane_dir, 0.0)) {
                params.tStep = -num_voxels / 2;
            } else {
                params.tStep = num_voxels / 2;
            }
            params.within_bounds = true;
            return params;
        }
        if (t_min_within_bounds && (t_min < t_max || isKnEqual(t, t_max))) {
            params.tStep = -1;
            params.tMax = t_min;
            params.within_bounds = true;
            return params;
        }
        if (KnLessThan(t, t_max) && KnLessThan(t_max, t_end) && (t_max < t_min || isKnEqual(t, t_min))) {
            params.tStep = 1;
            params.tMax = t_max;
            params.within_bounds = true;
            return params;
        }
    }
    params.tStep = 0;
    params.tMax = std::numeric_limits<double>::infinity();
    params.within_bounds = false;
    return params;
}

// Determines whether an angular hit occurs for the given ray. An angular hit is considered an intersection with
// the ray and an angular section. The angular sections live in the XY plane.
// Input:
//    ray: The given ray to check for intersection.
//    grid: The grid that the ray is intersecting with.
//    p*_angular_*: Points of intersection between the lines corresponding to angular voxels.
//    boundaries and the initial radial voxel of the ray.
//    t: The current time.
//    t_end: The maximum allowable time before the minimum of (sphere exit, grid exit, t_end).
//
// Returns: the corresponding angular hit parameters.
AngularHitParameters angularHit(const Ray &ray, const SphericalVoxelGrid &grid, double px_angular_one,
                                double px_angular_two, double py_angular_one, double py_angular_two,
                                double t, double t_end) noexcept {
    // Ray segment vector.
    const BoundVec3 p = ray.pointAtParameter(t);
    const BoundVec3 p_end = ray.pointAtParameter(t_end);
    const FreeVec3 v = p_end - p;

    // Calculate the voxel boundary vectors.
    const FreeVec3 p_one(px_angular_one, py_angular_one, 0.0);
    const FreeVec3 p_two(px_angular_two, py_angular_two, 0.0);
    const BoundVec3 u_min = grid.sphereCenter() - p_one;
    const BoundVec3 u_max = grid.sphereCenter() - p_two;
    const FreeVec3 w_min = p_one - FreeVec3(p);
    const FreeVec3 w_max = p_two - FreeVec3(p);
    const double perp_uv_min = u_min.x() * v.y() - u_min.y() * v.x();
    const double perp_uv_max = u_max.x() * v.y() - u_max.y() * v.x();
    const double perp_uw_min = u_min.x() * w_min.y() - u_min.y() * w_min.x();
    const double perp_uw_max = u_max.x() * w_max.y() - u_max.y() * w_max.x();
    const double perp_vw_min = v.x() * w_min.y() - v.y() * w_min.x();
    const double perp_vw_max = v.x() * w_max.y() - v.y() * w_max.x();
    const GenHitParameters params = generalizedPlaneHit(ray, perp_uv_min, perp_uv_max, perp_uw_min, perp_uw_max,
                                                        perp_vw_min, perp_vw_max, p, v, t, t_end, ray.direction().y(),
                                                        grid.numAngularVoxels());
    return {.tMaxTheta=params.tMax, .tStepTheta=params.tStep, .within_bounds=params.within_bounds};
}

// Determines whether an azimuthal hit occurs for the given ray. An azimuthal hit is considered an intersection with
// the ray and an azimuthal section. The azimuthal sections live in the XZ plane.
// Input:
//    ray: The given ray to check for intersection.
//    grid: The grid that the ray is intersecting with.
//    p*_azimuthal_*: Points of intersection between the lines corresponding to azimuthal voxels.
//    boundaries and the initial radial voxel of the ray.
//    t: The current time.
//    t_end: The maximum allowable time before the minimum of (sphere exit, grid exit, t_end).
//    v: The dot product between the ray's direction and the ray_sphere_vector.
//
// Returns: the corresponding azimuthal hit parameters.
AzimuthalHitParameters azimuthalHit(const Ray &ray, const SphericalVoxelGrid &grid, double px_azimuthal_one,
                                    double px_azimuthal_two, double pz_azimuthal_one, double pz_azimuthal_two,
                                    double t, double t_end) noexcept {
    // Ray segment vector.
    const BoundVec3 p = ray.pointAtParameter(t);
    const BoundVec3 p_end = ray.pointAtParameter(t_end);
    const FreeVec3 v = p_end - p;

    // Calculate the voxel boundary vectors.
    const FreeVec3 p_one(px_azimuthal_one, 0.0, pz_azimuthal_one);
    const FreeVec3 p_two(px_azimuthal_two, 0.0, pz_azimuthal_two);
    const BoundVec3 u_min = grid.sphereCenter() - p_one;
    const BoundVec3 u_max = grid.sphereCenter() - p_two;
    const FreeVec3 w_min = p_one - FreeVec3(p);
    const FreeVec3 w_max = p_two - FreeVec3(p);
    const double perp_uv_min = u_min.x() * v.z() - u_min.z() * v.x();
    const double perp_uv_max = u_max.x() * v.z() - u_max.z() * v.x();
    const double perp_uw_min = u_min.x() * w_min.z() - u_min.z() * w_min.x();
    const double perp_uw_max = u_max.x() * w_max.z() - u_max.z() * w_max.x();
    const double perp_vw_min = v.x() * w_min.z() - v.z() * w_min.x();
    const double perp_vw_max = v.x() * w_max.z() - v.z() * w_max.x();
    const GenHitParameters params = generalizedPlaneHit(ray, perp_uv_min, perp_uv_max, perp_uw_min, perp_uw_max,
                                                        perp_vw_min, perp_vw_max, p, v, t, t_end, ray.direction().z(),
                                                        grid.numAzimuthalVoxels());
    return {.tMaxPhi=params.tMax, .tStepPhi=params.tStep, .within_bounds=params.within_bounds};
}

// Calculates the voxel(s) with the minimal tMax for the next intersection. Since t is being updated
// with each interval of the algorithm, this must check the following cases:
// 1. tMaxTheta is the minimum.
// 2. tMaxR is the minimum.
// 3. tMaxPhi is the minimum.
// 4. tMaxTheta, tMaxR, tMaxPhi equal intersection.
// 5. tMaxR, tMaxPhi equal intersection.
// 6. tMaxR, tMaxTheta equal intersection.
// 7. tMaxPhi, tMaxTheta equal intersection.
//
// For each case, the following must hold:
// (1) t < tMax < t_end
// (2) If its a radial hit, the next step must be within bounds.
//
// In cases 1, 3 we also need to change the requirements slightly for when
// the ray intersects a single shell, but still crosses an angular or azimuthal boundary:
// (1) tMax must be within bounds
// (2) Either tMax is a strict minimum OR the next step is a radial exit.
inline VoxelIntersectionType minimumIntersection(const RadialHitParameters &rad_params,
                                                 const AngularHitParameters &ang_params,
                                                 const AzimuthalHitParameters &azi_params) noexcept {
    if (ang_params.within_bounds && ((ang_params.tMaxTheta < rad_params.tMaxR
                                      && rad_params.tMaxR < azi_params.tMaxPhi) || rad_params.exits_voxel_bounds)) {
        return VoxelIntersectionType::Angular;
    }
    if (rad_params.within_bounds && rad_params.tMaxR < ang_params.tMaxTheta
        && rad_params.tMaxR < azi_params.tMaxPhi) {
        return VoxelIntersectionType::Radial;
    }
    if (azi_params.within_bounds && ((azi_params.tMaxPhi < ang_params.tMaxTheta
                                      && azi_params.tMaxPhi < rad_params.tMaxR) || rad_params.exits_voxel_bounds)) {
        return VoxelIntersectionType::Azimuthal;
    }
    if (rad_params.within_bounds && isKnEqual(rad_params.tMaxR, ang_params.tMaxTheta)
        && isKnEqual(rad_params.tMaxR, azi_params.tMaxPhi)) {
        return VoxelIntersectionType::RadialAngularAzimuthal;
    }
    if (azi_params.within_bounds && isKnEqual(azi_params.tMaxPhi, ang_params.tMaxTheta)) {
        return VoxelIntersectionType::AngularAzimuthal;
    }
    if (rad_params.within_bounds && isKnEqual(ang_params.tMaxTheta, rad_params.tMaxR)) {
        return VoxelIntersectionType::RadialAngular;
    }
    if (rad_params.within_bounds && isKnEqual(rad_params.tMaxR, azi_params.tMaxPhi)) {
        return VoxelIntersectionType::RadialAzimuthal;
    }
    return VoxelIntersectionType::None;
}

// Create an array of values representing the points of intersection between the lines corresponding
// to angular voxel boundaries and the initial radial voxel of the ray in the XY plane. This is similar for
// azimuthal voxel boundaries, but in the XZ plane instead.
inline void initializeVoxelBoundarySegments(std::vector<double> &Px_angular, std::vector<double> &Py_angular,
                                            std::vector<double> &Px_azimuthal, std::vector<double> &Pz_azimuthal,
                                            const SphericalVoxelGrid &grid, double current_r) noexcept {
    double radians = 0;
    for (std::size_t j = 0; j < Px_angular.size(); ++j) {
        Px_angular[j] = current_r * std::cos(radians) + grid.sphereCenter().x();
        Py_angular[j] = current_r * std::sin(radians) + grid.sphereCenter().y();
        radians += grid.deltaTheta();
    }
    radians = 0;
    for (std::size_t n = 0; n < Px_azimuthal.size(); ++n) {
        Px_azimuthal[n] = current_r * std::cos(radians) + grid.sphereCenter().x();
        Pz_azimuthal[n] = current_r * std::sin(radians) + grid.sphereCenter().z();
        radians += grid.deltaPhi();
    }
}

// Each time the radius changes, we need to update the angular and azimuthal voxel boundary segments.
// This function updates Px_angular, Py_angular, Px_azimuthal, and Pz_azimuthal.
inline void updateVoxelBoundarySegments(std::vector<double> &Px_angular, std::vector<double> &Py_angular,
                                        std::vector<double> &Px_azimuthal, std::vector<double> &Pz_azimuthal,
                                        const SphericalVoxelGrid &grid, std::size_t current_voxel_ID_r) noexcept {
    const double new_r = grid.sphereMaxRadius() - grid.deltaRadius() * (current_voxel_ID_r - 1);
    for (std::size_t l = 0; l < Px_angular.size(); ++l) {
        const double old_angular_x = grid.sphereCenter().x() - Px_angular[l];
        const double old_angular_y = grid.sphereCenter().y() - Py_angular[l];
        const double new_r_over_plane_length = new_r / std::sqrt(old_angular_x * old_angular_x +
                                                                 old_angular_y * old_angular_y);
        Px_angular[l] = grid.sphereCenter().x() - new_r_over_plane_length * old_angular_x;
        Py_angular[l] = grid.sphereCenter().y() - new_r_over_plane_length * old_angular_y;
    }
    for (std::size_t m = 0; m < Px_azimuthal.size(); ++m) {
        const double old_azimuthal_x = grid.sphereCenter().x() - Px_azimuthal[m];
        const double old_azimuthal_z = grid.sphereCenter().z() - Pz_azimuthal[m];
        const double new_r_over_plane_length = new_r / std::sqrt(old_azimuthal_x * old_azimuthal_x +
                                                                 old_azimuthal_z * old_azimuthal_z);
        Px_azimuthal[m] = grid.sphereCenter().x() - new_r_over_plane_length * old_azimuthal_x;
        Pz_azimuthal[m] = grid.sphereCenter().z() - new_r_over_plane_length * old_azimuthal_z;
    }
}

// p1 will lie between two angular voxel boundaries iff the angle between it and the angular boundary intersection
// points along the circle of max radius is obtuse. Equality represents the case when the point lies on an angular
// boundary. This is similar for azimuthal boundaries. Since both cases use points in a plane (XY for angular, XZ
// for azimuthal), this can be generalized to a single function.
inline int calculateVoxelID(const std::vector<double> plane1, const std::vector<double> plane2,
                            double p1, double p2) noexcept {
    std::size_t i = 0;
    while (i < plane1.size() - 1) {
        const double px_diff = plane1[i] - plane1[i + 1];
        const double py_diff = plane2[i] - plane2[i + 1];
        const double px_p1_diff = plane1[i] - p1;
        const double py_p1_diff = plane2[i] - p2;
        const double n_px_p1_diff = plane1[i + 1] - p1;
        const double n_py_p1_diff = plane2[i + 1] - p2;
        const double d1 = (px_p1_diff * px_p1_diff) + (py_p1_diff * py_p1_diff);
        const double d2 = (n_px_p1_diff * n_px_p1_diff) + (n_py_p1_diff * n_py_p1_diff);
        const double d3 = (px_diff * px_diff) + (py_diff * py_diff);
        if (d1 + d2 <= d3) { return i; }
        ++i;
    }
    return -1;
}

std::vector<SphericalVoxel> sphericalCoordinateVoxelTraversal(const Ray &ray, const SphericalVoxelGrid &grid,
                                                              double t_begin, double t_end) noexcept {
    std::vector<SphericalVoxel> voxels;
    voxels.reserve(grid.numRadialVoxels() + grid.numAngularVoxels() + grid.numAzimuthalVoxels());
    /* INITIALIZATION PHASE */

    // Determine ray location at t_begin.
    const BoundVec3 point_at_t_begin = ray.pointAtParameter(t_begin);
    const FreeVec3 ray_sphere_vector = grid.sphereCenter() - point_at_t_begin;
    double entry_radius = grid.deltaRadius(); // The radius at which the ray first enters the sphere.

    const double ray_sphere_vector_dot = ray_sphere_vector.dot(ray_sphere_vector);
    while (ray_sphere_vector_dot > (entry_radius * entry_radius) && entry_radius < grid.sphereMaxRadius()) {
        entry_radius += grid.deltaRadius();
    }

    // Find the intersection times for the ray and the radial shell containing the parameter point at t_begin.
    // This will determine if the ray intersects the sphere.
    const double v = ray_sphere_vector.dot(ray.unitDirection().to_free());
    const double discriminant = (entry_radius * entry_radius) - (ray_sphere_vector_dot - v * v);

    if (discriminant <= 0.0) { return voxels; }
    const double d = std::sqrt(discriminant);

    // Calculate the time of entrance and exit of the ray.
    // Need to use a non-zero direction to determine this.
    const double t_begin_t1 = ray.timeOfIntersectionAt(v - d);
    const double t_begin_t2 = ray.timeOfIntersectionAt(v + d);

    if ((t_begin_t1 < t_begin && t_begin_t2 < t_begin) || isKnEqual(t_begin_t1, t_begin_t2)) {
        // Case 1: No intersection.
        // Case 2: Tangent hit.
        return voxels;
    }
    int current_voxel_ID_r = 1 + (grid.sphereMaxRadius() - entry_radius) * grid.invDeltaRadius();

    std::vector<double> Px_angular(grid.numAngularVoxels() + 1);
    std::vector<double> Py_angular(grid.numAngularVoxels() + 1);
    std::vector<double> Px_azimuthal(grid.numAzimuthalVoxels() + 1);
    std::vector<double> Pz_azimuthal(grid.numAzimuthalVoxels() + 1);
    initializeVoxelBoundarySegments(Px_angular, Py_angular, Px_azimuthal, Pz_azimuthal, grid, entry_radius);

    double a, b, c;
    const bool ray_origin_is_outside_grid = isKnEqual(entry_radius, grid.sphereMaxRadius());
    if (isKnEqual(ray.origin(), grid.sphereCenter())) {
        // If the ray starts at the sphere's center, we need to perturb slightly along
        // the path to determine the correct angular and azimuthal voxel.
        const double perturbed_t = 0.1;
        a = grid.sphereCenter().x() - (point_at_t_begin.x() + ray.direction().x() * perturbed_t);
        b = grid.sphereCenter().y() - (point_at_t_begin.y() + ray.direction().y() * perturbed_t);
        c = grid.sphereCenter().z() - (point_at_t_begin.z() + ray.direction().z() * perturbed_t);
    } else if (ray_origin_is_outside_grid) {
        const BoundVec3 pa = ray.pointAtParameter(t_begin_t1);
        a = grid.sphereCenter().x() - pa.x();
        b = grid.sphereCenter().y() - pa.y();
        c = grid.sphereCenter().z() - pa.z();
    } else {
        a = grid.sphereCenter().x() - ray.origin().x();
        b = grid.sphereCenter().y() - ray.origin().y();
        c = grid.sphereCenter().z() - ray.origin().z();
    }

    double p_x, p_y, p_z;
    const double ang_plane_length = std::sqrt(a * a + b * b);
    if (isKnEqual(ang_plane_length, 0.0)) {
        p_x = grid.sphereCenter().x() + entry_radius;
        p_y = grid.sphereCenter().y();
    } else {
        const double r_over_ang_plane_length = entry_radius / ang_plane_length;
        p_x = grid.sphereCenter().x() - a * r_over_ang_plane_length;
        p_y = grid.sphereCenter().y() - b * r_over_ang_plane_length;
    }
    int current_voxel_ID_theta = calculateVoxelID(Px_angular, Py_angular, p_x, p_y);

    const double azi_plane_length = std::sqrt(a * a + c * c);
    if (isKnEqual(azi_plane_length, 0.0)) {
        p_x = grid.sphereCenter().x() + entry_radius;
        p_z = grid.sphereCenter().z();
    } else {
        const double r_over_azi_plane_length = entry_radius / azi_plane_length;
        p_x = grid.sphereCenter().x() - a * r_over_azi_plane_length;
        p_z = grid.sphereCenter().z() - c * r_over_azi_plane_length;
    }
    int current_voxel_ID_phi = calculateVoxelID(Px_azimuthal, Pz_azimuthal, p_x, p_z);

    voxels.push_back({.radial_voxel=current_voxel_ID_r,
                             .angular_voxel=current_voxel_ID_theta,
                             .azimuthal_voxel=current_voxel_ID_phi});

    // Find the maximum time the ray will be in the grid.
    const double max_discriminant = grid.sphereMaxRadius() * grid.sphereMaxRadius() - (ray_sphere_vector_dot - v * v);
    const double max_d = std::sqrt(max_discriminant);
    const double t_grid_exit = std::max(ray.timeOfIntersectionAt(v - max_d), ray.timeOfIntersectionAt(v + max_d));

    // Find the correct time to begin the traversal phase.
    double t;
    if (ray_origin_is_outside_grid) {
        t = ray.timeOfIntersectionAt(Vec3(p_x, p_y, p_z));
    } else { t = t_begin; }

    /* TRAVERSAL PHASE */
    bool previous_transition_flag = false;
    bool radius_has_stepped = false;
    t_end = std::min(t_grid_exit, t_end);

    while (t < t_end) {
        const auto radial_params = radialHit(ray, grid, current_voxel_ID_r, ray_sphere_vector_dot,
                                             t, t_end, v, previous_transition_flag);
        previous_transition_flag = radial_params.previous_transition_flag;
        const auto angular_params = angularHit(ray, grid, Px_angular[current_voxel_ID_theta],
                                               Px_angular[current_voxel_ID_theta + 1],
                                               Py_angular[current_voxel_ID_theta],
                                               Py_angular[current_voxel_ID_theta + 1], t, t_end);
        const auto azimuthal_params = azimuthalHit(ray, grid, Px_azimuthal[current_voxel_ID_phi],
                                                   Px_azimuthal[current_voxel_ID_phi + 1],
                                                   Pz_azimuthal[current_voxel_ID_phi],
                                                   Pz_azimuthal[current_voxel_ID_phi + 1], t, t_end);
        const auto voxel_intersection = minimumIntersection(radial_params, angular_params, azimuthal_params);
        switch (voxel_intersection) {
            case Radial: {
                t = radial_params.tMaxR;
                current_voxel_ID_r += radial_params.tStepR;
                radius_has_stepped = radial_params.tStepR != 0;
                break;
            }
            case Angular: {
                t = angular_params.tMaxTheta;
                current_voxel_ID_theta = (current_voxel_ID_theta + angular_params.tStepTheta) % grid.numAngularVoxels();
                break;
            }
            case Azimuthal: {
                t = azimuthal_params.tMaxPhi;
                current_voxel_ID_phi = (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                break;
            }
            case RadialAngular: {
                t = radial_params.tMaxR;
                current_voxel_ID_r += radial_params.tStepR;
                current_voxel_ID_theta = (current_voxel_ID_theta + angular_params.tStepTheta) % grid.numAngularVoxels();
                radius_has_stepped = radial_params.tStepR != 0;
                break;
            }
            case RadialAzimuthal: {
                t = radial_params.tMaxR;
                current_voxel_ID_r += radial_params.tStepR;
                current_voxel_ID_phi = (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                radius_has_stepped = radial_params.tStepR != 0;
                break;
            }
            case AngularAzimuthal: {
                t = angular_params.tMaxTheta;
                current_voxel_ID_theta = (current_voxel_ID_theta + angular_params.tStepTheta) % grid.numAngularVoxels();
                current_voxel_ID_phi = (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                break;
            }
            case RadialAngularAzimuthal: {
                t = radial_params.tMaxR;
                current_voxel_ID_r += radial_params.tStepR;
                current_voxel_ID_theta = (current_voxel_ID_theta + angular_params.tStepTheta) % grid.numAngularVoxels();
                current_voxel_ID_phi = (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                radius_has_stepped = radial_params.tStepR != 0;
                break;
            }
            case None: {
                return voxels;
            }
        }
        if (radius_has_stepped) {
            updateVoxelBoundarySegments(Px_angular, Py_angular, Px_azimuthal, Pz_azimuthal, grid, current_voxel_ID_r);
            radius_has_stepped = false;
        }
        voxels.push_back({.radial_voxel=current_voxel_ID_r,
                                 .angular_voxel=current_voxel_ID_theta,
                                 .azimuthal_voxel=current_voxel_ID_phi});
    }
    return voxels;
}

std::vector<SphericalVoxel> sphericalCoordinateVoxelTraversalCy(double *ray_origin, double *ray_direction,
                                                                double *min_bound, double *max_bound,
                                                                std::size_t num_radial_voxels,
                                                                std::size_t num_angular_voxels,
                                                                std::size_t num_azimuthal_voxels, double *sphere_center,
                                                                double sphere_max_radius, double t_begin,
                                                                double t_end) noexcept {
    const Ray ray(BoundVec3(ray_origin[0], ray_origin[1], ray_origin[2]),
                  FreeVec3(ray_direction[0], ray_direction[1], ray_direction[2]));
    const SphericalVoxelGrid grid(BoundVec3(min_bound[0], min_bound[1], min_bound[2]),
                                  BoundVec3(max_bound[0], max_bound[1], max_bound[2]),
                                  num_radial_voxels, num_angular_voxels, num_azimuthal_voxels,
                                  BoundVec3(sphere_center[0], sphere_center[1], sphere_center[2]), sphere_max_radius);
    return sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
}