#include "spherical_volume_rendering_util.h"
#include <vector>
#include <array>
#include <algorithm>
#include <limits>

namespace svr {
    // An array of step values used to avoid branch prediction in radialHit().
    constexpr std::array<int, 3> step{0, -1, 1};

    // Epsilons used for floating point comparisons in Knuth's algorithm.
    constexpr double ABS_EPSILON = 1e-12;
    constexpr double REL_EPSILON = 1e-8;

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

    // The necessary information to calculate a radial hit. ray_sphere_vector_dot and v are calculated in the
    // initialization phase, so unnecessary to re-calculate again for radial hit. intersection_times and times_gt_t
    // are used to determine the type of radial hit. Upon initialization, the previous transition flag is set to false.
    struct RadialHitData {
        inline RadialHitData(double t_v, double t_ray_sphere_vector_dot) :
        v(t_v), ray_sphere_vector_dot(t_ray_sphere_vector_dot),
        previous_transition_flag(false) { times_gt_t.reserve(4); }

        // Pre-calculated data to be used when calculating a radial hit.
        double v, ray_sphere_vector_dot;
        // Pre-initialized structures to be used when calculating a radial hit.
        std::array<double, 4> intersection_times;
        std::vector<double> times_gt_t;

        // The current state of the previous_transition_flag. This is saved here so that it can be passed
        // into the radial hit function with each traversal.
        bool previous_transition_flag;
    };

    // Pre-calculated information for the generalized plane hit function, which generalizes azimuthal and angular
    // hits. Since the ray segment is dependent solely on time, this is unnecessary to calculate twice for each
    // plane hit function. Here, ray_segment is the difference between P2 and P1.
    // The collinear times are the two times possible for t, dependent on if the ray is collinear to the given
    // voxel boundary.
    struct RaySegment {
        inline RaySegment(const Ray &ray, double t_end) : P2(ray.pointAtParameter(t_end)) {}

        // Updates the point P1 with the new time traversal time t. Similarly, updates the
        // segment denoted by P2 - P1.
        inline void updateAtTime(double t, const Ray &ray) noexcept {
            P1 = ray.pointAtParameter(t);
            ray_segment = P2 - P1;
        }

        // The begin point of the ray segment.
        BoundVec3 P1 = BoundVec3(0.0, 0.0, 0.0);

        // The end point of the ray segment.
        BoundVec3 P2 = BoundVec3(0.0, 0.0, 0.0);

        // P2 - P1.
        FreeVec3 ray_segment = FreeVec3(0.0, 0.0, 0.0);
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

    // A point will lie between two angular voxel boundaries iff the angle between it and the angular boundary
    // intersection points along the circle of max radius is obtuse. Equality represents the case when the point lies
    // on an angular boundary. This is similar for azimuthal boundaries. Since both cases use points in a plane
    // (XY for angular, XZ for azimuthal), this can be generalized to a single function. Since angular and azimuthal
    // voxel boundaries range from [0, N], returns -1 in the case where the point does not lie within the boundaries.
    inline int calculateVoxelID(const std::vector<svr::LineSegment> &plane, double p1, double p2) noexcept {
        std::size_t i = 0;
        while (i < plane.size() - 1) {
            const double px_diff = plane[i].P1 - plane[i + 1].P1;
            const double py_diff = plane[i].P2 - plane[i + 1].P2;
            const double px_p1_diff = plane[i].P1 - p1;
            const double py_p1_diff = plane[i].P2 - p2;
            const double n_px_p1_diff = plane[i + 1].P1 - p1;
            const double n_py_p1_diff = plane[i + 1].P2 - p2;
            const double d1d2 = (px_p1_diff * px_p1_diff) + (py_p1_diff * py_p1_diff) +
                                (n_px_p1_diff * n_px_p1_diff) + (n_py_p1_diff * n_py_p1_diff);
            const double d3 = (px_diff * px_diff) + (py_diff * py_diff);
            if (KnLessThan(d1d2, d3) || isKnEqual(d1d2, d3)) { return i; }
            ++i;
        }
        return -1;
    }

    // Determines whether a radial hit occurs for the given ray. A radial hit is considered an intersection with
    // the ray and a radial section. The struct RadialHitData is used to provide already initialized data structures,
    // as well as avoiding unnecessary duplicate calculations that have already been done in the initialization phase.
    // This follows closely the mathematics presented in:
    // http://cas.xav.free.fr/Graphics%20Gems%204%20-%20Paul%20S.%20Heckbert.pdf
    RadialHitParameters radialHit(const Ray &ray, const svr::SphericalVoxelGrid &grid, RadialHitData &rdata,
                                  int current_voxel_ID_r, double t, double t_end) noexcept {
        const double current_radius = grid.sphereMaxRadius() - grid.deltaRadius() * (current_voxel_ID_r - 1);
        double r_a = std::max(current_radius - grid.deltaRadius(), grid.deltaRadius());
        double r_b;
        if (!rdata.previous_transition_flag) {
            // To find the next radius, we need to check the previous_transition_flag:
            // In the case that the ray has sequential hits with equal radii, e.g.
            // the innermost radial disc, this ensures that the proper radii are being checked.
            r_b = std::min(current_radius + grid.deltaRadius(), grid.sphereMaxRadius());
        } else { r_b = std::min(current_radius, grid.sphereMaxRadius()); }
        // Find the intersection times for the ray and the previous and next radial discs.
        const double ray_sphere_dot_minus_v_squared = rdata.ray_sphere_vector_dot - rdata.v * rdata.v;
        double discriminant_a = r_a * r_a - ray_sphere_dot_minus_v_squared;
        if (discriminant_a < 0.0) {
            r_a += grid.deltaRadius();
            discriminant_a = r_a * r_a - ray_sphere_dot_minus_v_squared;
        }
        const double d_a = std::sqrt(discriminant_a);

        rdata.intersection_times[0] = ray.timeOfIntersectionAt(rdata.v - d_a);
        rdata.intersection_times[1] = ray.timeOfIntersectionAt(rdata.v + d_a);
        rdata.intersection_times[2] = 0.0;
        rdata.intersection_times[3] = 0.0;

        const double discriminant_b = r_b * r_b - ray_sphere_dot_minus_v_squared;
        if (discriminant_b >= 0.0) {
            const double d_b = std::sqrt(discriminant_b);
            rdata.intersection_times[2] = ray.timeOfIntersectionAt(rdata.v - d_b);
            rdata.intersection_times[3] = ray.timeOfIntersectionAt(rdata.v + d_b);
        }
        rdata.times_gt_t.clear();
        std::copy_if(rdata.intersection_times.cbegin(), rdata.intersection_times.cend(),
                     std::back_inserter(rdata.times_gt_t), [t](double i) { return i > t; });

        if (rdata.times_gt_t.empty()) {
            // No intersection.
            return {.tMaxR=std::numeric_limits<double>::max(),
                    .tStepR=0,
                    .previous_transition_flag=false,
                    .within_bounds=false};
        }
        const double intersection_time = rdata.times_gt_t[0];
        const double r_new = (ray.pointAtParameter(intersection_time) - grid.sphereCenter()).length();
        const bool is_radial_transition = isKnEqual(r_new, current_radius);
        const bool is_not_tangential_hit = !(rdata.times_gt_t.size() >= 2 &&
                                             isKnEqual(rdata.intersection_times[0], rdata.intersection_times[1]));
        return {.tMaxR=intersection_time,
                .tStepR=step[1 * is_not_tangential_hit + (is_not_tangential_hit &&          // { 0, -1, 1 }
                             !is_radial_transition && KnLessThan(r_new, current_radius))],
                .previous_transition_flag=is_radial_transition,
                .within_bounds=KnLessThan(t, intersection_time) && KnLessThan(intersection_time, t_end)
        };
    }

    // A generalized version of the latter half of the angular and azimuthal hit parameters. Since the only difference
    // is the 2-d plane for which they exist in, this portion can be generalized to a single function call.
    // The variables that are generalized take the form of *_plane_*, such as ray_plane_direction. If this called in
    // AngularHit(), ray_plane_direction == ray.direction.y(). The calculations presented below follow closely the
    // works of [Foley et al, 1996], [O'Rourke, 1998].
    // Reference: http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
    GenHitParameters
    generalizedPlaneHit(const svr::SphericalVoxelGrid &grid, const Ray &ray, double perp_uv_min, double perp_uv_max,
                        double perp_uw_min, double perp_uw_max, double perp_vw_min, double perp_vw_max,
                        const RaySegment &rs_data, const std::array<double, 2> &collinear_times, double t, double t_end,
                        double ray_plane_direction, double sphere_plane_center,
                        const std::vector<svr::LineSegment> &P_max, int current_voxel_ID) noexcept {
        const bool is_parallel_min = isKnEqual(perp_uv_min, 0.0);
        const bool is_collinear_min = is_parallel_min && isKnEqual(perp_uw_min, 0.0) && isKnEqual(perp_vw_min, 0.0);
        const bool is_parallel_max = isKnEqual(perp_uv_max, 0.0);
        const bool is_collinear_max = is_parallel_max && isKnEqual(perp_uw_max, 0.0) && isKnEqual(perp_vw_max, 0.0);
        double a, b;
        double t_min = collinear_times[is_collinear_min];
        bool is_intersect_min = false;
        if (!is_parallel_min) {
            const double inv_perp_uv_min = 1.0 / perp_uv_min;
            a = perp_vw_min * inv_perp_uv_min;
            b = perp_uw_min * inv_perp_uv_min;
            if (!((KnLessThan(a, 0.0) || KnLessThan(1.0, a)) || KnLessThan(b, 0.0) || KnLessThan(1.0, b))) {
                is_intersect_min = true;
                t_min = ray.timeOfIntersectionAt(Vec3(rs_data.P1.x() + rs_data.ray_segment.x() * b,
                                                      rs_data.P1.y() + rs_data.ray_segment.y() * b,
                                                      rs_data.P1.z() + rs_data.ray_segment.z() * b));
            }
        }
        double t_max = collinear_times[is_collinear_max];
        bool is_intersect_max = false;
        if (!is_parallel_max) {
            const double inv_perp_uv_max = 1.0 / perp_uv_max;
            a = perp_vw_max * inv_perp_uv_max;
            b = perp_uw_max * inv_perp_uv_max;
            if (!((KnLessThan(a, 0.0) || KnLessThan(1.0, a)) || KnLessThan(b, 0.0) || KnLessThan(1.0, b))) {
                is_intersect_max = true;
                t_max = ray.timeOfIntersectionAt(Vec3(rs_data.P1.x() + rs_data.ray_segment.x() * b,
                                                      rs_data.P1.y() + rs_data.ray_segment.y() * b,
                                                      rs_data.P1.z() + rs_data.ray_segment.z() * b));
            }
        }
        GenHitParameters params;
        if (is_intersect_max && !is_intersect_min && !is_collinear_min && KnLessThan(t_max, t_end) &&
            KnLessThan(t, t_max)) {
            params.tStep = 1;
            params.tMax = t_max;
            params.within_bounds = true;
            return params;
        }
        if (is_intersect_min && !is_intersect_max && !is_collinear_max && KnLessThan(t_min, t_end) &&
            KnLessThan(t, t_min)) {
            params.tStep = -1;
            params.tMax = t_min;
            params.within_bounds = true;
            return params;
        }
        if ((is_intersect_min && is_intersect_max) ||
            (is_intersect_min && is_collinear_max) ||
            (is_intersect_max && is_collinear_min)) {
            const bool t_min_within_bounds = KnLessThan(t, t_min) && KnLessThan(t_min, t_end);
            if (t_min_within_bounds && isKnEqual(t_min, t_max)) {
                params.tMax = t_max;
                const double perturbed_t = 0.1;
                a = -ray.direction().x() * perturbed_t;
                b = -ray_plane_direction * perturbed_t;
                const double max_radius_over_plane_length = grid.sphereMaxRadius() / std::sqrt(a * a + b * b);
                const double p1 = grid.sphereCenter().x() - max_radius_over_plane_length * a;
                const double p2 = sphere_plane_center - max_radius_over_plane_length * b;
                const int next_step = std::abs(current_voxel_ID - calculateVoxelID(P_max, p1, p2));

                params.tStep = (KnLessThan(ray_plane_direction, 0.0) || KnLessThan(ray.direction().x(), 0.0)) ?
                               next_step : -next_step;
                params.within_bounds = true;
                return params;
            }
            if (t_min_within_bounds && (KnLessThan(t_min, t_max) || isKnEqual(t, t_max))) {
                params.tStep = -1;
                params.tMax = t_min;
                params.within_bounds = true;
                return params;
            }
            const bool t_max_within_bounds = KnLessThan(t, t_max) && KnLessThan(t_max, t_end);
            if (t_max_within_bounds && (KnLessThan(t_max, t_min) || isKnEqual(t, t_min))) {
                params.tStep = 1;
                params.tMax = t_max;
                params.within_bounds = true;
                return params;
            }
        }
        params.tStep = 0;
        params.tMax = std::numeric_limits<double>::max();
        params.within_bounds = false;
        return params;
    }

    // Determines whether an angular hit occurs for the given ray. An angular hit is considered an intersection with
    // the ray and an angular section. The angular sections live in the XY plane.
    AngularHitParameters angularHit(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                    const RaySegment &rs_data, const std::array<double, 2> &collinear_times,
                                    int current_voxel_ID_theta, double t, double t_end) noexcept {
        // Calculate the voxel boundary vectors.
        const FreeVec3 p_one(grid.pMaxAngular(current_voxel_ID_theta).P1,
                             grid.pMaxAngular(current_voxel_ID_theta).P2, 0.0);
        const FreeVec3 p_two(grid.pMaxAngular(current_voxel_ID_theta + 1).P1,
                             grid.pMaxAngular(current_voxel_ID_theta + 1).P2, 0.0);
        const BoundVec3 u_min = grid.sphereCenter() - p_one;
        const BoundVec3 u_max = grid.sphereCenter() - p_two;
        const FreeVec3 w_min = p_one - FreeVec3(rs_data.P1);
        const FreeVec3 w_max = p_two - FreeVec3(rs_data.P1);
        const double perp_uv_min = u_min.x() * rs_data.ray_segment.y() - u_min.y() * rs_data.ray_segment.x();
        const double perp_uv_max = u_max.x() * rs_data.ray_segment.y() - u_max.y() * rs_data.ray_segment.x();
        const double perp_uw_min = u_min.x() * w_min.y() - u_min.y() * w_min.x();
        const double perp_uw_max = u_max.x() * w_max.y() - u_max.y() * w_max.x();
        const double perp_vw_min = rs_data.ray_segment.x() * w_min.y() - rs_data.ray_segment.y() * w_min.x();
        const double perp_vw_max = rs_data.ray_segment.x() * w_max.y() - rs_data.ray_segment.y() * w_max.x();
        const GenHitParameters params = generalizedPlaneHit(grid, ray, perp_uv_min, perp_uv_max, perp_uw_min,
                                                            perp_uw_max, perp_vw_min, perp_vw_max, rs_data,
                                                            collinear_times, t, t_end, ray.direction().y(),
                                                            grid.sphereCenter().y(), grid.pMaxAngular(),
                                                            current_voxel_ID_theta);
        return {.tMaxTheta=params.tMax, .tStepTheta=params.tStep, .within_bounds=params.within_bounds};
    }

    // Determines whether an azimuthal hit occurs for the given ray. An azimuthal hit is considered an intersection with
    // the ray and an azimuthal section. The azimuthal sections live in the XZ plane.
    AzimuthalHitParameters azimuthalHit(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                        const RaySegment &rs_data, const std::array<double, 2> &collinear_times,
                                        int current_voxel_ID_phi, double t, double t_end) noexcept {
        // Calculate the voxel boundary vectors.
        const FreeVec3 p_one(grid.pMaxAzimuthal(current_voxel_ID_phi).P1, 0.0,
                             grid.pMaxAzimuthal(current_voxel_ID_phi).P2);
        const FreeVec3 p_two(grid.pMaxAzimuthal(current_voxel_ID_phi + 1).P1, 0.0,
                             grid.pMaxAzimuthal(current_voxel_ID_phi + 1).P2);
        const BoundVec3 u_min = grid.sphereCenter() - p_one;
        const BoundVec3 u_max = grid.sphereCenter() - p_two;
        const FreeVec3 w_min = p_one - FreeVec3(rs_data.P1);
        const FreeVec3 w_max = p_two - FreeVec3(rs_data.P1);
        const double perp_uv_min = u_min.x() * rs_data.ray_segment.z() - u_min.z() * rs_data.ray_segment.x();
        const double perp_uv_max = u_max.x() * rs_data.ray_segment.z() - u_max.z() * rs_data.ray_segment.x();
        const double perp_uw_min = u_min.x() * w_min.z() - u_min.z() * w_min.x();
        const double perp_uw_max = u_max.x() * w_max.z() - u_max.z() * w_max.x();
        const double perp_vw_min = rs_data.ray_segment.x() * w_min.z() - rs_data.ray_segment.z() * w_min.x();
        const double perp_vw_max = rs_data.ray_segment.x() * w_max.z() - rs_data.ray_segment.z() * w_max.x();
        const GenHitParameters params = generalizedPlaneHit(grid, ray, perp_uv_min, perp_uv_max, perp_uw_min,
                                                            perp_uw_max, perp_vw_min, perp_vw_max, rs_data,
                                                            collinear_times, t, t_end, ray.direction().z(),
                                                            grid.sphereCenter().z(), grid.pMaxAzimuthal(),
                                                            current_voxel_ID_phi);
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
    inline VoxelIntersectionType minimumIntersection(const RadialHitParameters &rad_params,
                                                     const AngularHitParameters &ang_params,
                                                     const AzimuthalHitParameters &azi_params) noexcept {
        if (rad_params.within_bounds && KnLessThan(rad_params.tMaxR, ang_params.tMaxTheta)
            && KnLessThan(rad_params.tMaxR, azi_params.tMaxPhi)) {
            return VoxelIntersectionType::Radial;
        }
        if (ang_params.within_bounds && KnLessThan(ang_params.tMaxTheta, rad_params.tMaxR)
            && KnLessThan(ang_params.tMaxTheta, azi_params.tMaxPhi)) {
            return VoxelIntersectionType::Angular;
        }
        if (azi_params.within_bounds && KnLessThan(azi_params.tMaxPhi, ang_params.tMaxTheta)
            && KnLessThan(azi_params.tMaxPhi, rad_params.tMaxR)) {
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

    // Initialize an array of values representing the points of intersection between the lines corresponding
    // to voxel boundaries and a given radial voxel in the XY plane and XZ plane. Here, P_* represents
    // these points with a given radius 'current_radius'. The case where the number of angular voxels is
    // equal to the number of azimuthal voxels is also checked to reduce the number of trigonometric
    // and floating point calculations.
    inline void initializeVoxelBoundarySegments(std::vector<svr::LineSegment> &P_angular,
                                                std::vector<svr::LineSegment> &P_azimuthal,
                                                const svr::SphericalVoxelGrid &grid, double current_radius) noexcept {
        if (grid.numAngularVoxels() == grid.numAzimuthalVoxels()) {
            std::transform(grid.angularTrigValues().cbegin(), grid.angularTrigValues().cend(),
                           P_angular.begin(), P_azimuthal.begin(),
                           [current_radius, &grid](const TrigonometricValues &tv,
                                                   LineSegment &angular_LS) -> LineSegment {
                               const double px_value = current_radius * tv.cosine + grid.sphereCenter().x();
                               const double current_radius_times_sin = current_radius * tv.sine;
                               angular_LS = {.P1=px_value, .P2=current_radius_times_sin + grid.sphereCenter().y()};
                               return {.P1=px_value, .P2=current_radius_times_sin + grid.sphereCenter().z()};
                           });
            return;
        }
        std::transform(grid.angularTrigValues().cbegin(), grid.angularTrigValues().cend(), P_angular.begin(),
                       [current_radius, &grid](const TrigonometricValues &tv) -> LineSegment {
                           return {.P1=current_radius * tv.cosine + grid.sphereCenter().x(),
                                   .P2=current_radius * tv.sine + grid.sphereCenter().y()};
                       });
        std::transform(grid.azimuthalTrigValues().cbegin(), grid.azimuthalTrigValues().cend(), P_azimuthal.begin(),
                       [current_radius, &grid](const TrigonometricValues &tv) -> LineSegment {
                           return {.P1=current_radius * tv.cosine + grid.sphereCenter().x(),
                                   .P2=current_radius * tv.sine + grid.sphereCenter().z()};
                       });
    }

    std::vector<svr::SphericalVoxel> sphericalCoordinateVoxelTraversal(const Ray &ray,
                                                                       const svr::SphericalVoxelGrid &grid,
                                                                       double t_begin, double t_end) noexcept {
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

        if (discriminant <= 0.0) { return {}; }
        const double d = std::sqrt(discriminant);

        // Calculate the time of entrance and exit of the ray.
        // Need to use a non-zero direction to determine this.
        const double t_begin_t1 = ray.timeOfIntersectionAt(v - d);
        const double t_begin_t2 = ray.timeOfIntersectionAt(v + d);

        if ((t_begin_t1 < t_begin && t_begin_t2 < t_begin) || isKnEqual(t_begin_t1, t_begin_t2)) {
            // Case 1: No intersection.
            // Case 2: Tangent hit.
            return {};
        }
        int current_voxel_ID_r = 1 + (grid.sphereMaxRadius() - entry_radius) * grid.invDeltaRadius();
        std::vector<svr::LineSegment> P_angular(grid.numAngularVoxels() + 1);
        std::vector<svr::LineSegment> P_azimuthal(grid.numAzimuthalVoxels() + 1);
        initializeVoxelBoundarySegments(P_angular, P_azimuthal, grid, entry_radius);

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
        int current_voxel_ID_theta = calculateVoxelID(P_angular, p_x, p_y);

        const double azi_plane_length = std::sqrt(a * a + c * c);
        if (isKnEqual(azi_plane_length, 0.0)) {
            p_x = grid.sphereCenter().x() + entry_radius;
            p_z = grid.sphereCenter().z();
        } else {
            const double r_over_azi_plane_length = entry_radius / azi_plane_length;
            p_x = grid.sphereCenter().x() - a * r_over_azi_plane_length;
            p_z = grid.sphereCenter().z() - c * r_over_azi_plane_length;
        }
        int current_voxel_ID_phi = calculateVoxelID(P_azimuthal, p_x, p_z);

        std::vector<svr::SphericalVoxel> voxels;
        voxels.reserve(grid.numRadialVoxels() + grid.numAngularVoxels() + grid.numAzimuthalVoxels());
        voxels.push_back({.radial_voxel=current_voxel_ID_r,
                          .angular_voxel=current_voxel_ID_theta,
                          .azimuthal_voxel=current_voxel_ID_phi});

        // Find the maximum time the ray will be in the grid.
        const double max_discriminant =
                grid.sphereMaxRadius() * grid.sphereMaxRadius() - (ray_sphere_vector_dot - v * v);
        const double max_d = std::sqrt(max_discriminant);
        const double t_grid_exit = std::max(ray.timeOfIntersectionAt(v - max_d), ray.timeOfIntersectionAt(v + max_d));
        // Find the correct time to begin the traversal phase.
        double t = ray_origin_is_outside_grid ? ray.timeOfIntersectionAt(Vec3(p_x, p_y, p_z)) : t_begin;
        t_end = std::min(t_grid_exit, t_end);

        // Initialize the time in case of collinear min or collinear max for generalized plane hits.
        std::array<double, 2> collinear_times = {0.0, ray.timeOfIntersectionAt(grid.sphereCenter())};

        RadialHitData radial_hit_data(v, ray_sphere_vector_dot);
        RaySegment ray_segment(ray, t_end);

        while (true) {
            const auto radial_params = radialHit(ray, grid, radial_hit_data, current_voxel_ID_r, t, t_end);
            radial_hit_data.previous_transition_flag = radial_params.previous_transition_flag;
            ray_segment.updateAtTime(t, ray);
            const auto angular_params = angularHit(ray, grid, ray_segment, collinear_times,
                                                   current_voxel_ID_theta, t, t_end);
            const auto azimuthal_params = azimuthalHit(ray, grid, ray_segment, collinear_times,
                                                       current_voxel_ID_phi, t, t_end);
            const auto voxel_intersection = minimumIntersection(radial_params, angular_params, azimuthal_params);
            switch (voxel_intersection) {
                case Radial: {
                    t = radial_params.tMaxR;
                    current_voxel_ID_r += radial_params.tStepR;
                    break;
                }
                case Angular: {
                    t = angular_params.tMaxTheta;
                    current_voxel_ID_theta =
                            (current_voxel_ID_theta + angular_params.tStepTheta) % grid.numAngularVoxels();
                    break;
                }
                case Azimuthal: {
                    t = azimuthal_params.tMaxPhi;
                    current_voxel_ID_phi =
                            (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                    break;
                }
                case RadialAngular: {
                    t = radial_params.tMaxR;
                    current_voxel_ID_r += radial_params.tStepR;
                    current_voxel_ID_theta =
                            (current_voxel_ID_theta + angular_params.tStepTheta) % grid.numAngularVoxels();
                    break;
                }
                case RadialAzimuthal: {
                    t = radial_params.tMaxR;
                    current_voxel_ID_r += radial_params.tStepR;
                    current_voxel_ID_phi =
                            (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                    break;
                }
                case AngularAzimuthal: {
                    t = angular_params.tMaxTheta;
                    current_voxel_ID_theta =
                            (current_voxel_ID_theta + angular_params.tStepTheta) % grid.numAngularVoxels();
                    current_voxel_ID_phi =
                            (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                    break;
                }
                case RadialAngularAzimuthal: {
                    t = radial_params.tMaxR;
                    current_voxel_ID_r += radial_params.tStepR;
                    current_voxel_ID_theta =
                            (current_voxel_ID_theta + angular_params.tStepTheta) % grid.numAngularVoxels();
                    current_voxel_ID_phi =
                            (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                    break;
                }
                case None: {
                    return voxels;
                }
            }
            voxels.push_back({.radial_voxel=current_voxel_ID_r,
                              .angular_voxel=current_voxel_ID_theta,
                              .azimuthal_voxel=current_voxel_ID_phi});
        }
        return voxels;
    }

    std::vector<svr::SphericalVoxel> sphericalCoordinateVoxelTraversalCy(double *ray_origin, double *ray_direction,
                                                                         double *min_bound, double *max_bound,
                                                                         std::size_t num_radial_voxels,
                                                                         std::size_t num_angular_voxels,
                                                                         std::size_t num_azimuthal_voxels,
                                                                         double *sphere_center,
                                                                         double sphere_max_radius, double t_begin,
                                                                         double t_end) noexcept {
        const Ray ray(BoundVec3(ray_origin[0], ray_origin[1], ray_origin[2]),
                      FreeVec3(ray_direction[0], ray_direction[1], ray_direction[2]));
        const svr::SphericalVoxelGrid grid(BoundVec3(min_bound[0], min_bound[1], min_bound[2]),
                                           BoundVec3(max_bound[0], max_bound[1], max_bound[2]),
                                           num_radial_voxels, num_angular_voxels, num_azimuthal_voxels,
                                           BoundVec3(sphere_center[0], sphere_center[1], sphere_center[2]),
                                           sphere_max_radius);
        return svr::sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
    }

} // namespace svr