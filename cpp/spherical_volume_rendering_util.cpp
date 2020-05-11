#include "spherical_volume_rendering_util.h"
#include <vector>
#include <array>
#include <algorithm>
#include <limits>

namespace svr {
    // An array of step values used to avoid branch prediction in radialHit().
    constexpr std::array<int, 3> STEP{0, -1, 1};

    // Represents an invalid time. At no point should the time be a negative number.
    constexpr double INVALID_TIME = -1.0;

    // Epsilon used for floating point comparisons in Knuth's algorithm.
    constexpr double ABS_EPSILON = 1e-12;
    constexpr double REL_EPSILON = 1e-8;

    // The type corresponding to the voxel(s) with the minimum tMax value for a given traversal.
    enum VoxelIntersectionType {
        None = 0,
        Radial = 1,
        Polar = 2,
        Azimuthal = 3,
        RadialPolar = 4,
        RadialAzimuthal = 5,
        PolarAzimuthal = 6,
        RadialPolarAzimuthal = 7
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

    // The parameters returned by polarHit().
    struct PolarHitParameters {
        // The time at which a hit occurs for the ray at the next point of intersection with a polar section.
        // This is always calculated from the ray origin.
        double tMaxTheta;

        // The voxel traversal value of a polar step: 0, +1, -1. This is added to the current polar voxel.
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

    // A generalized set of angular hit parameters.
    struct AngularHitParameters {
        double tMax;
        int tStep;
        bool within_bounds;
    };

    // A lot of calculations are conducted in the initialization phase already. To mitigate these from occurring again,
    // This structure allows the calculations to be saved for use during radialHit().
    struct RadialHitData {
    public:
        inline RadialHitData(double v, double rsvd_minus_v_squared) :
                v_(v), rsvd_minus_v_squared_(rsvd_minus_v_squared) {}

        inline double v() const noexcept { return v_; }

        inline double rsvdMinusVSquared() const noexcept { return rsvd_minus_v_squared_;}

        inline bool transitionFlag() const noexcept { return transition_flag_; }

        inline void updateTransitionFlag(bool b) noexcept { transition_flag_ = b; }

    private:
        // Pre-calculated data to be used when calculating a radial hit.
        const double v_, rsvd_minus_v_squared_;

        // The current state of the previous_transition_flag. It marks if we've made a transition
        // in sign for radial steps, i.e. + => - or - => +.
        bool transition_flag_;
    };

    // Pre-calculated information for the generalized plane hit function, which generalizes azimuthal and polar
    // hits. Since the ray segment is dependent solely on time, this is unnecessary to calculate twice for each
    // plane hit function. Here, ray_segment is the difference between P2 and P1.
    // The collinear times are the two times possible for t, dependent on if the ray is collinear to the given
    // voxel boundary.
    struct RaySegment {
    public:
        inline RaySegment(const Ray *ray, double t_end) : ray_(ray), P2_(ray->pointAtParameter(t_end)) {}

        // Updates the point P1 with the new time traversal time t. Similarly, updates the
        // segment denoted by P2 - P1.
        inline void updateAtTime(double t) noexcept {
            P1_ = ray_->pointAtParameter(t);
            ray_segment_ = P2_ - P1_;
        }

        // Calculates the updated ray segment intersection point given an intersect parameter.
        // More information on this use case can be found at:
        // http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
        inline double intersectionTimeAt(double intersect_param) const noexcept {
            return (P1_[NZDI_] + ray_segment_[NZDI_] * intersect_param - ray_->origin()[NZDI_])
                   * ray_->invDirection()[NZDI_];
        }

        inline const BoundVec3 &P1() const noexcept { return P1_; }

        inline const BoundVec3 &P2() const noexcept { return P2_; }

        inline const FreeVec3 &raySegment() const noexcept { return ray_segment_; }

    private:
        // The associated ray for which the segment (P1, P2) refers to.
        const Ray *ray_;

        // The non-zero direction index of the ray.
        const DirectionIndex NZDI_ = ray_->NonZeroDirectionIndex();

        // The begin point of the ray segment.
        BoundVec3 P1_;

        // The end point of the ray segment.
        const BoundVec3 P2_;

        // P2 - P1.
        FreeVec3 ray_segment_;
    };

    // Determines equality between two floating point numbers using an absolute epsilon.
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

    // Checks to see if a is strictly less than b with an absolute epsilon.
    inline bool lessThan(double a, double b) noexcept {
        return a < b && !isEqual(a, b);
    }

    // Maintains a count until BinaryPredicate p is met for two consecutive elements.
    template<class ForwardIt, class BinaryPredicate>
    inline std::size_t index_adjacent_find_until(ForwardIt first, ForwardIt last, BinaryPredicate p) noexcept {
        std::size_t count = 0;
        ForwardIt next = first;
        ++next;
        for (; next != last; ++next, ++first) {
            if (p(*first, *next)) return count;
            count++;
        }
        return count;
    }

    // A point will lie between two polar voxel boundaries iff the angle between it and the polar boundary
    // intersection points along the circle of max radius is obtuse. Equality represents the case when the point lies
    // on an polar boundary. This is similar for azimuthal boundaries. Since both cases use points in a plane
    // (XY for polar, XZ for azimuthal), this can be generalized to a single function.
    inline int calculateVoxelIDFromPoints(const std::vector<LineSegment> &angular_max, const double p1, double p2) noexcept {
        return index_adjacent_find_until(angular_max.cbegin(), angular_max.cend(),
                                         [p1, p2](const LineSegment &LS1, const LineSegment &LS2) -> bool {
                                             const double X_diff = LS1.P1 - LS2.P1;
                                             const double Y_diff = LS1.P2 - LS2.P2;
                                             const double X_p1_diff = LS1.P1 - p1;
                                             const double X_p2_diff = LS1.P2 - p2;
                                             const double Y_p1_diff = LS2.P1 - p1;
                                             const double Y_p2_diff = LS2.P2 - p2;
                                             const double d1d2 = (X_p1_diff * X_p1_diff) + (X_p2_diff * X_p2_diff) +
                                                                 (Y_p1_diff * Y_p1_diff) + (Y_p2_diff * Y_p2_diff);
                                             const double d3 = (X_diff * X_diff) + (Y_diff * Y_diff);
                                             return d1d2 < d3 || isEqual(d1d2, d3);
                                         });
    }


    // Initializes an angular voxel ID. For polar initializations=, *_2 represents the y-plane. For azimuthal
    // initialization, it represents the z-plane. If the squared euclidean distance of the ray_sphere vector in the
    // given plane is zero, the voxel ID is set to 0.
    inline int initializeAngularVoxelID(const SphericalVoxelGrid &grid, const FreeVec3 &ray_sphere,
                                        const std::vector<LineSegment> &angular_max,
                                        double ray_sphere_2, double grid_sphere_2, double entry_radius) noexcept {
        const double SED = ray_sphere.x() * ray_sphere.x() + ray_sphere_2 * ray_sphere_2;
        if (isEqual(SED, 0.0)) { return 0; }
        const double r = entry_radius / std::sqrt(SED);
        const double p1 = grid.sphereCenter().x() - ray_sphere.x() * r;
        const double p2 = grid_sphere_2 - ray_sphere_2 * r;
        return calculateVoxelIDFromPoints(angular_max, p1, p2);
    }

    // Determines whether a radial hit occurs for the given ray. A radial hit is considered an intersection with
    // the ray and a radial section. The struct RadialHitData is used to provide already initialized data structures,
    // as well as avoiding unnecessary duplicate calculations that have already been done in the initialization phase.
    // This follows closely the mathematics presented in:
    // http://cas.xav.free.fr/Graphics%20Gems%204%20-%20Paul%20S.%20Heckbert.pdf
    inline RadialHitParameters radialHit(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                         const RadialHitData &rh_data, int current_voxel_ID_r,
                                         double t, double t_end) noexcept {
        const std::size_t voxel_idx = current_voxel_ID_r - 1;
        const double current_radius = grid.deltaRadii(voxel_idx);

        // Find the intersection times for the ray and the previous radial disc.
        const std::size_t previous_idx = std::min(voxel_idx + 1, grid.numRadialVoxels() - 1);
        const double r_a = grid.deltaRadiiSquared(previous_idx -
                                                  (grid.deltaRadiiSquared(previous_idx) < rh_data.rsvdMinusVSquared()));
        const double d_a = std::sqrt(r_a - rh_data.rsvdMinusVSquared());

        std::array<double, 4> intersection_times;
        intersection_times[0] = ray.timeOfIntersectionAt(rh_data.v() - d_a);
        intersection_times[1] = ray.timeOfIntersectionAt(rh_data.v() + d_a);

        // To find the next radius, we need to check the previous_transition_flag:
        // In the case that the ray has sequential hits with equal radii, e.g.
        // the innermost radial disc, this ensures that the proper radii are being checked.
        const double transition_radii[] = {grid.deltaRadiiSquared( std::min(voxel_idx - 1, std::size_t{0}) ),
                                           grid.deltaRadiiSquared( voxel_idx ) };
        const double r_b = transition_radii[rh_data.transitionFlag()];
        if (r_b >= rh_data.rsvdMinusVSquared()) {
            const double d_b = std::sqrt(r_b - rh_data.rsvdMinusVSquared());
            intersection_times[2] = ray.timeOfIntersectionAt(rh_data.v() - d_b);
            intersection_times[3] = ray.timeOfIntersectionAt(rh_data.v() + d_b);
        }
        const auto intersection_time_it = std::find_if(intersection_times.cbegin(), intersection_times.cend(),
                                                       [t](double i)->double{ return i > t; });
        if (intersection_time_it == intersection_times.cend()) {
            return {.tMaxR=std::numeric_limits<double>::max(), .tStepR=0,
                    .previous_transition_flag=false, .within_bounds=false };
        }
        const double intersection_time = *intersection_time_it;
        const double r_new = (ray.pointAtParameter(intersection_time) - grid.sphereCenter()).length();
        const bool is_radial_transition = isEqual(r_new, current_radius);
        const bool is_tangential_hit = isEqual(intersection_times[0], intersection_times[1]);
        return {.tMaxR=intersection_time,
                .tStepR=STEP[1 * !is_tangential_hit +
                             (!is_tangential_hit && !is_radial_transition
                              && lessThan(r_new, current_radius))],
                .previous_transition_flag=is_radial_transition,
                .within_bounds=lessThan(t, intersection_time) && lessThan(intersection_time, t_end) };
    }

    // A generalized version of the latter half of the polar and azimuthal hit parameters. Since the only difference
    // is the 2-d plane for which they exist in, this portion can be generalized to a single function call.
    // The calculations presented below follow closely the works of [Foley et al, 1996], [O'Rourke, 1998].
    // Reference: http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()

    AngularHitParameters
    angularHit(const svr::SphericalVoxelGrid &grid, const Ray &ray, double perp_uv_min, double perp_uv_max,
               double perp_uw_min, double perp_uw_max, double perp_vw_min, double perp_vw_max,
               const RaySegment &RS, const std::array<double, 2> &collinear_times, double t, double t_end,
               double ray_direction_2, double sphere_center_2,
               const std::vector<svr::LineSegment> &P_max, int current_voxel_ID) noexcept {
        const bool is_parallel_min = isEqual(perp_uv_min, 0.0);
        const bool is_collinear_min = is_parallel_min && isEqual(perp_uw_min, 0.0) && isEqual(perp_vw_min, 0.0);
        const bool is_parallel_max = isEqual(perp_uv_max, 0.0);
        const bool is_collinear_max = is_parallel_max && isEqual(perp_uw_max, 0.0) && isEqual(perp_vw_max, 0.0);
        double a, b;
        double t_min = collinear_times[is_collinear_min];
        bool is_intersect_min = false;
        if (!is_parallel_min) {
            const double inv_perp_uv_min = 1.0 / perp_uv_min;
            a = perp_vw_min * inv_perp_uv_min;
            b = perp_uw_min * inv_perp_uv_min;
            if (!((lessThan(a, 0.0) || lessThan(1.0, a)) || lessThan(b, 0.0) || lessThan(1.0, b))) {
                is_intersect_min = true;
                t_min = RS.intersectionTimeAt(b);
            }
        }
        double t_max = collinear_times[is_collinear_max];
        bool is_intersect_max = false;
        if (!is_parallel_max) {
            const double inv_perp_uv_max = 1.0 / perp_uv_max;
            a = perp_vw_max * inv_perp_uv_max;
            b = perp_uw_max * inv_perp_uv_max;
            if (!((lessThan(a, 0.0) || lessThan(1.0, a)) || lessThan(b, 0.0) || lessThan(1.0, b))) {
                is_intersect_max = true;
                t_max = RS.intersectionTimeAt(b);
            }
        }

        const bool t_max_within_bounds = lessThan(t, t_max) && lessThan(t_max, t_end);
        const bool t_min_within_bounds = lessThan(t, t_min) && lessThan(t_min, t_end);
        if (!t_max_within_bounds && !t_min_within_bounds) {
            return { .tMax = std::numeric_limits<double>::max(), .tStep = 0, .within_bounds = false };
        }


        if (t_max_within_bounds && is_intersect_max && !is_intersect_min && !is_collinear_min) {
            return { .tMax = t_max, .tStep = 1, .within_bounds = true };
        }
        if (t_min_within_bounds && is_intersect_min && !is_intersect_max && !is_collinear_max) {
            return { .tMax = t_min, .tStep = -1, .within_bounds = true };
        }
        if ((is_intersect_min && is_intersect_max) ||
            (is_intersect_min && is_collinear_max) ||
            (is_intersect_max && is_collinear_min)) {
            if (t_min_within_bounds && isEqual(t_min, t_max)) {
                const double perturbed_t = 0.1;
                a = -ray.direction().x() * perturbed_t;
                b = -ray_direction_2 * perturbed_t;
                const double max_radius_over_plane_length = grid.sphereMaxRadius() / std::sqrt(a * a + b * b);
                const double p1 = grid.sphereCenter().x() - max_radius_over_plane_length * a;
                const double p2 = sphere_center_2 - max_radius_over_plane_length * b;
                const int next_step = std::abs(current_voxel_ID - calculateVoxelIDFromPoints(P_max, p1, p2));
                return {.tMax = t_max,
                        .tStep = (lessThan(ray_direction_2, 0.0) || lessThan(ray.direction().x(), 0.0)) ?
                                 next_step : -next_step,
                        .within_bounds = true
                };
            }
            if (t_min_within_bounds && (lessThan(t_min, t_max) || isEqual(t, t_max))) {
                return { .tMax = t_min, .tStep = -1, .within_bounds = true };
            }
            if (t_max_within_bounds && (lessThan(t_max, t_min) || isEqual(t, t_min))) {
                return { .tMax = t_max, .tStep = 1, .within_bounds = true };
            }
        }
        return { .tMax = std::numeric_limits<double>::max(), .tStep = 0, .within_bounds = false };
    }

    // Determines whether a polar hit occurs for the given ray. A polar hit is considered an intersection with
    // the ray and a polar section. The polar sections live in the XY plane.
    inline PolarHitParameters polarHit(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                       const RaySegment &RS, const std::array<double, 2> &collinear_times,
                                       int current_voxel_ID_theta, double t, double t_end) noexcept {
        // Calculate the voxel boundary vectors.
        const FreeVec3 p_one(grid.pMaxPolar(current_voxel_ID_theta).P1,
                             grid.pMaxPolar(current_voxel_ID_theta).P2, 0.0);
        const FreeVec3 p_two(grid.pMaxPolar(current_voxel_ID_theta + 1).P1,
                             grid.pMaxPolar(current_voxel_ID_theta + 1).P2, 0.0);
        const BoundVec3 u_min = grid.sphereCenter() - p_one;
        const BoundVec3 u_max = grid.sphereCenter() - p_two;
        const FreeVec3 w_min = p_one - FreeVec3(RS.P1());
        const FreeVec3 w_max = p_two - FreeVec3(RS.P1());
        const double perp_uv_min = u_min.x() * RS.raySegment().y() - u_min.y() * RS.raySegment().x();
        const double perp_uv_max = u_max.x() * RS.raySegment().y() - u_max.y() * RS.raySegment().x();
        const double perp_uw_min = u_min.x() * w_min.y() - u_min.y() * w_min.x();
        const double perp_uw_max = u_max.x() * w_max.y() - u_max.y() * w_max.x();
        const double perp_vw_min = RS.raySegment().x() * w_min.y() - RS.raySegment().y() * w_min.x();
        const double perp_vw_max = RS.raySegment().x() * w_max.y() - RS.raySegment().y() * w_max.x();
        const AngularHitParameters params = angularHit(grid, ray, perp_uv_min, perp_uv_max, perp_uw_min,
                                                       perp_uw_max, perp_vw_min, perp_vw_max, RS,
                                                       collinear_times, t, t_end, ray.direction().y(),
                                                       grid.sphereCenter().y(), grid.pMaxPolar(),
                                                       current_voxel_ID_theta);
        return {.tMaxTheta=params.tMax, .tStepTheta=params.tStep, .within_bounds=params.within_bounds};
    }

    // Determines whether an azimuthal hit occurs for the given ray. An azimuthal hit is
    // considered an intersection with the ray and an azimuthal section.
    // The azimuthal sections live in the XZ plane.
    inline AzimuthalHitParameters azimuthalHit(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                               const RaySegment &RS, const std::array<double, 2> &collinear_times,
                                               int current_voxel_ID_phi, double t, double t_end) noexcept {
        // Calculate the voxel boundary vectors.
        const FreeVec3 p_one(grid.pMaxAzimuthal(current_voxel_ID_phi).P1, 0.0,
                             grid.pMaxAzimuthal(current_voxel_ID_phi).P2);
        const FreeVec3 p_two(grid.pMaxAzimuthal(current_voxel_ID_phi + 1).P1, 0.0,
                             grid.pMaxAzimuthal(current_voxel_ID_phi + 1).P2);
        const BoundVec3 u_min = grid.sphereCenter() - p_one;
        const BoundVec3 u_max = grid.sphereCenter() - p_two;
        const FreeVec3 w_min = p_one - FreeVec3(RS.P1());
        const FreeVec3 w_max = p_two - FreeVec3(RS.P1());
        const double perp_uv_min = u_min.x() * RS.raySegment().z() - u_min.z() * RS.raySegment().x();
        const double perp_uv_max = u_max.x() * RS.raySegment().z() - u_max.z() * RS.raySegment().x();
        const double perp_uw_min = u_min.x() * w_min.z() - u_min.z() * w_min.x();
        const double perp_uw_max = u_max.x() * w_max.z() - u_max.z() * w_max.x();
        const double perp_vw_min = RS.raySegment().x() * w_min.z() - RS.raySegment().z() * w_min.x();
        const double perp_vw_max = RS.raySegment().x() * w_max.z() - RS.raySegment().z() * w_max.x();
        const AngularHitParameters params = angularHit(grid, ray, perp_uv_min, perp_uv_max, perp_uw_min,
                                                       perp_uw_max, perp_vw_min, perp_vw_max, RS,
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
                                                     const PolarHitParameters &ang_params,
                                                     const AzimuthalHitParameters &azi_params) noexcept {
        if (rad_params.within_bounds && lessThan(rad_params.tMaxR, ang_params.tMaxTheta)
            && lessThan(rad_params.tMaxR, azi_params.tMaxPhi)) {
            return VoxelIntersectionType::Radial;
        }
        if (ang_params.within_bounds && lessThan(ang_params.tMaxTheta, rad_params.tMaxR)
            && lessThan(ang_params.tMaxTheta, azi_params.tMaxPhi)) {
            return VoxelIntersectionType::Polar;
        }
        if (azi_params.within_bounds && lessThan(azi_params.tMaxPhi, ang_params.tMaxTheta)
            && lessThan(azi_params.tMaxPhi, rad_params.tMaxR)) {
            return VoxelIntersectionType::Azimuthal;
        }
        if (rad_params.within_bounds && isEqual(rad_params.tMaxR, ang_params.tMaxTheta)
            && isEqual(rad_params.tMaxR, azi_params.tMaxPhi)) {
            return VoxelIntersectionType::RadialPolarAzimuthal;
        }
        if (azi_params.within_bounds && isEqual(azi_params.tMaxPhi, ang_params.tMaxTheta)) {
            return VoxelIntersectionType::PolarAzimuthal;
        }
        if (rad_params.within_bounds && isEqual(ang_params.tMaxTheta, rad_params.tMaxR)) {
            return VoxelIntersectionType::RadialPolar;
        }
        if (rad_params.within_bounds && isEqual(rad_params.tMaxR, azi_params.tMaxPhi)) {
            return VoxelIntersectionType::RadialAzimuthal;
        }
        return VoxelIntersectionType::None;
    }

    // Initialize an array of values representing the points of intersection between the lines corresponding
    // to voxel boundaries and a given radial voxel in the XY plane and XZ plane. Here, P_* represents
    // these points with a given radius 'current_radius'. The case where the number of polar voxels is
    // equal to the number of azimuthal voxels is also checked to reduce the number of trigonometric
    // and floating point calculations.
    inline void initializeVoxelBoundarySegments(std::vector<svr::LineSegment> &P_polar,
                                                std::vector<svr::LineSegment> &P_azimuthal,
                                                bool ray_origin_is_outside_grid,
                                                const svr::SphericalVoxelGrid &grid, double current_radius) noexcept {
        if (ray_origin_is_outside_grid) {
            P_polar = grid.pMaxPolar();
            P_azimuthal = grid.pMaxAzimuthal();
            return;
        }

        if (grid.numPolarVoxels() == grid.numAzimuthalVoxels()) {
            std::transform(grid.polarTrigValues().cbegin(), grid.polarTrigValues().cend(),
                           P_polar.begin(), P_azimuthal.begin(),
                           [current_radius, &grid](const TrigonometricValues &tv,
                                                   LineSegment &polar_LS) -> LineSegment {
                               const double px_value = current_radius * tv.cosine + grid.sphereCenter().x();
                               const double current_radius_times_sin = current_radius * tv.sine;
                               polar_LS = {.P1=px_value, .P2=current_radius_times_sin + grid.sphereCenter().y()};
                               return {.P1=px_value, .P2=current_radius_times_sin + grid.sphereCenter().z()};
                           });
            return;
        }
        std::transform(grid.polarTrigValues().cbegin(), grid.polarTrigValues().cend(), P_polar.begin(),
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

    std::vector<svr::SphericalVoxel> walkSphericalVolume(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                                         const svr::SphereBound &min_bound,
                                                         const svr::SphereBound &max_bound,
                                                         double t_begin, double t_end) noexcept {
        const FreeVec3 rsv = grid.sphereCenter() - ray.pointAtParameter(t_begin);   // Ray Sphere Vector.
        const FreeVec3 rsv_tz = (t_begin == 0.0) ? rsv : grid.sphereCenter() - ray.pointAtParameter(0.0);

        const double rsvd_begin = rsv.dot(rsv);
        std::size_t idx = grid.numRadialVoxels();
        const auto it = std::find_if(grid.deltaRadiiSquared().crbegin() + 1, grid.deltaRadiiSquared().crend(),
                                     [rsvd_begin, &idx](double dR_squared)-> bool {
                                         --idx;
                                         return rsvd_begin <= dR_squared;
                                     });
        const bool ray_origin_is_outside_grid = (it == grid.deltaRadiiSquared().crend());
        const double max_radius_squared = grid.deltaRadiiSquared()[0];
        const double entry_radius_squared = ray_origin_is_outside_grid ? max_radius_squared : *it;
        const double entry_radius = grid.deltaRadii()[idx];

        // Find the intersection times for the ray and the radial shell containing the parameter point at t_begin.
        // This will determine if the ray intersects the sphere.
        const double rsvd = rsv_tz.dot(rsv_tz); // Ray Sphere Vector Dot product at time zero.
        const double v = rsv_tz.dot(ray.unitDirection().to_free());
        const double rsvd_minus_v_squared = rsvd - v * v;

        if (entry_radius_squared <= rsvd_minus_v_squared) { return {}; }
        const double d = std::sqrt(entry_radius_squared - rsvd_minus_v_squared);

        // Calculate the time of entrance and exit of the ray.
        const double t_entrance = ray.timeOfIntersectionAt(v - d);
        const double t_exit = ray.timeOfIntersectionAt(v + d);

        if ((t_entrance < t_begin && t_exit < t_begin) || isEqual(t_entrance, t_exit)) { return {}; }
        int current_voxel_ID_r = idx + 1;
        const int min_radial_ID = std::max(min_bound.radial * grid.invDeltaRadius(), 1.0);
        const int max_radial_ID = max_bound.radial * grid.invDeltaRadius();
        if (current_voxel_ID_r < min_radial_ID) { return {}; }

        std::vector<svr::LineSegment> P_polar(grid.numPolarVoxels() + 1);
        std::vector<svr::LineSegment> P_azimuthal(grid.numAzimuthalVoxels() + 1);
        initializeVoxelBoundarySegments(P_polar, P_azimuthal, ray_origin_is_outside_grid, grid, entry_radius);

        const FreeVec3 ray_sphere = ray_origin_is_outside_grid                                  ?
                                    grid.sphereCenter() - ray.pointAtParameter(t_entrance)      :
                                    isEqual(rsv, Vec3(0.0, 0.0, 0.0))                           ?
                                    grid.sphereCenter() - ray.pointAtParameter(t_begin + 0.1)   :   rsv;

        int current_voxel_ID_theta = initializeAngularVoxelID(grid, ray_sphere, P_polar, ray_sphere.y(),
                                                              grid.sphereCenter().y(), entry_radius);
        const int min_polar_ID = min_bound.polar * grid.invDeltaTheta();
        const int max_polar_ID = max_bound.polar * grid.invDeltaTheta();
        if (current_voxel_ID_theta < min_polar_ID || current_voxel_ID_theta > max_polar_ID - 1) { return {}; }

        int current_voxel_ID_phi= initializeAngularVoxelID(grid, ray_sphere, P_azimuthal, ray_sphere.z(),
                                                           grid.sphereCenter().z(), entry_radius);
        const int min_azimuthal_ID = min_bound.azimuthal * grid.invDeltaPhi();
        const int max_azimuthal_ID = max_bound.azimuthal * grid.invDeltaPhi();
        if (current_voxel_ID_phi < min_azimuthal_ID || current_voxel_ID_phi > max_azimuthal_ID - 1) { return {}; }

        std::vector<svr::SphericalVoxel> voxels;
        voxels.reserve(grid.numRadialVoxels() + grid.numPolarVoxels() + grid.numAzimuthalVoxels());
        voxels.push_back({.radial_voxel=current_voxel_ID_r,
                          .polar_voxel=current_voxel_ID_theta,
                          .azimuthal_voxel=current_voxel_ID_phi});

        double t;
        if (ray_origin_is_outside_grid) {
            t = t_entrance;
            t_end = std::min(t_end, t_exit);
        } else {
            t = t_begin;
            const double max_d = std::sqrt(max_radius_squared - rsvd_minus_v_squared);
            t_end = std::min(t_end, std::max(ray.timeOfIntersectionAt(v - max_d), ray.timeOfIntersectionAt(v + max_d)));
        }

        // Initialize the time in case of collinear min or collinear max for generalized plane hits.
        const std::array<double, 2> collinear_times = {INVALID_TIME, ray.timeOfIntersectionAt(grid.sphereCenter())};

        RadialHitData radial_hit_data(v, rsvd_minus_v_squared);
        RaySegment ray_segment(&ray, t_end);
        while (current_voxel_ID_r <= max_radial_ID &&
               current_voxel_ID_theta < max_polar_ID &&
               current_voxel_ID_phi < max_azimuthal_ID) {
            const auto radial_params = radialHit(ray, grid, radial_hit_data, current_voxel_ID_r, t, t_end);
            radial_hit_data.updateTransitionFlag(radial_params.previous_transition_flag);
            ray_segment.updateAtTime(t);
            const auto polar_params = polarHit(ray, grid, ray_segment, collinear_times,
                                               current_voxel_ID_theta, t, t_end);
            const auto azimuthal_params = azimuthalHit(ray, grid, ray_segment, collinear_times,
                                                       current_voxel_ID_phi, t, t_end);
            const auto voxel_intersection = minimumIntersection(radial_params, polar_params, azimuthal_params);
            switch (voxel_intersection) {
                case Radial: {
                    t = radial_params.tMaxR;
                    current_voxel_ID_r += radial_params.tStepR;
                    break;
                }
                case Polar: {
                    t = polar_params.tMaxTheta;
                    current_voxel_ID_theta =
                            (current_voxel_ID_theta + polar_params.tStepTheta) % grid.numPolarVoxels();
                    break;
                }
                case Azimuthal: {
                    t = azimuthal_params.tMaxPhi;
                    current_voxel_ID_phi =
                            (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                    break;
                }
                case RadialPolar: {
                    t = radial_params.tMaxR;
                    current_voxel_ID_r += radial_params.tStepR;
                    current_voxel_ID_theta =
                            (current_voxel_ID_theta + polar_params.tStepTheta) % grid.numPolarVoxels();
                    break;
                }
                case RadialAzimuthal: {
                    t = radial_params.tMaxR;
                    current_voxel_ID_r += radial_params.tStepR;
                    current_voxel_ID_phi =
                            (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                    break;
                }
                case PolarAzimuthal: {
                    t = polar_params.tMaxTheta;
                    current_voxel_ID_theta =
                            (current_voxel_ID_theta + polar_params.tStepTheta) % grid.numPolarVoxels();
                    current_voxel_ID_phi =
                            (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                    break;
                }
                case RadialPolarAzimuthal: {
                    t = radial_params.tMaxR;
                    current_voxel_ID_r += radial_params.tStepR;
                    current_voxel_ID_theta =
                            (current_voxel_ID_theta + polar_params.tStepTheta) % grid.numPolarVoxels();
                    current_voxel_ID_phi =
                            (current_voxel_ID_phi + azimuthal_params.tStepPhi) % grid.numAzimuthalVoxels();
                    break;
                }
                case None: {
                    return voxels;
                }
            }
            voxels.push_back({.radial_voxel=current_voxel_ID_r,
                              .polar_voxel=current_voxel_ID_theta,
                              .azimuthal_voxel=current_voxel_ID_phi});
        }
        return voxels;
    }

    std::vector<svr::SphericalVoxel> walkSphericalVolume(double *ray_origin, double *ray_direction,
                                                         std::size_t num_radial_voxels,
                                                         std::size_t num_polar_voxels,
                                                         std::size_t num_azimuthal_voxels,
                                                         double *sphere_center, double sphere_max_radius,
                                                         double *min_bound, double* max_bound,
                                                         double t_begin, double t_end) noexcept {
        return svr::walkSphericalVolume(Ray(BoundVec3(ray_origin[0], ray_origin[1], ray_origin[2]),
                                            FreeVec3(ray_direction[0], ray_direction[1], ray_direction[2])),
                                        svr::SphericalVoxelGrid(num_radial_voxels,
                                                                num_polar_voxels,
                                                                num_azimuthal_voxels,
                                                                BoundVec3(sphere_center[0], sphere_center[1],
                                                                          sphere_center[2]), sphere_max_radius),
                                        svr::SphereBound{.radial=min_bound[0],
                                                         .polar=min_bound[1],
                                                         .azimuthal=min_bound[2]},
                                        svr::SphereBound{.radial=max_bound[0],
                                                         .polar=max_bound[1],
                                                         .azimuthal=max_bound[2]},
                                        t_begin, t_end);
    }

} // namespace svr