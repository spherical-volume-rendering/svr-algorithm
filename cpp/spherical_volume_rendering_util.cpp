#include "spherical_volume_rendering_util.h"
#include "floating_point_comparison_util.h"
#include <vector>
#include <array>
#include <algorithm>
#include <limits>

namespace svr {
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
    struct HitParameters {
        // The time at which a hit occurs for the ray at the next point of intersection with a section.
        double tMax;

        // The voxel traversal value of a radial step: 0, +1, -1. This is added to the current voxel.
        int tStep;

        // Determines whether the current hit is within time bounds (t, t_end).
        bool within_bounds;
    };

    // Metadata for the radial hit parameters. There are two main components of the metadata: the transition flag and
    // the previous radial voxel. The transition flag is necessary to determine when tStepR's sign changes, and
    // the previous radial voxel ensures we do not hit the same voxel twice in the case of a tangential hit.
    struct RadialHitMetadata {
    public:
        inline bool radialStepTransitionHasOccurred() const noexcept { return transition_flag_; }

        inline void isRadialStepTransition(bool b) noexcept { transition_flag_ = b; }

        inline int previousRadialVoxel() const noexcept { return previous_radial_voxel_; }

        inline void updatePreviousRadialVoxel(int radial_voxel) noexcept { previous_radial_voxel_ = radial_voxel; }

    private:
        // Determine whether the current voxel traversal was a ' radial step transition'.
        // This is necessary to determine when the radial steps should go from negative to positive or vice-versa,
        // since the radial voxels go from 1..N..1, where N is the number of radial sections.
        bool transition_flag_ = false;

        // The previous radial voxel ID. This is necessary to maintain in the case of a tangential hit. If the current
        // radial voxel is equal to the previous radial voxel (i.e. a tangential hit occurred) and the only
        // voxel intersection is Radial, then we don't want another intersection to occur with the same voxel.
        int previous_radial_voxel_ = -1;
    };

    // Pre-calculated information for the generalized angular hit function, which generalizes azimuthal and polar
    // hits. Since the ray segment is dependent solely on time, this is unnecessary to calculate twice for each
    // plane hit function. Here, ray_segment is the difference between P2 and P1.
    struct RaySegment {
    public:
        inline RaySegment(double t_end, const Ray &ray) : P2_(ray.pointAtParameter(t_end)),
                                                          NZDI_(ray.NonZeroDirectionIndex()) {}

        // Updates the point P1 with the new time traversal time t. Similarly, updates the
        // segment denoted by P2 - P1.
        inline void updateAtTime(double t, const Ray &ray) noexcept {
            P1_ = ray.pointAtParameter(t);
            ray_segment_ = P2_ - P1_;
        }

        // Calculates the updated ray segment intersection point given an intersect parameter.
        // More information on the use case can be found at:
        // http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
        inline double intersectionTimeAt(double intersect_parameter, const Ray &ray) const noexcept {
            return (P1_[NZDI_] + ray_segment_[NZDI_] * intersect_parameter - ray.origin()[NZDI_])
                   * ray.invDirection()[NZDI_];
        }

        inline const BoundVec3 &P1() const noexcept { return P1_; }

        inline const BoundVec3 &P2() const noexcept { return P2_; }

        inline const FreeVec3 &vector() const noexcept { return ray_segment_; }

    private:
        // The end point of the ray segment.
        const BoundVec3 P2_;

        // The non-zero direction index of the ray.
        const DirectionIndex NZDI_;

        // The begin point of the ray segment.
        BoundVec3 P1_;

        // The free vector represented by P2 - P1.
        FreeVec3 ray_segment_;
    };

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
    inline int calculateAngularVoxelIDFromPoints(const std::vector<LineSegment> &angular_max,
                                                 const double p1, double p2) noexcept {
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
                                             return d1d2 < d3 || svr::isEqual(d1d2, d3);
                                         });
    }


    // Initializes an angular voxel ID. For polar initializations=, *_2 represents the y-plane. For azimuthal
    // initialization, it represents the z-plane. If the number of sections is 1 or the squared euclidean distance
    // of the ray_sphere vector in the given plane is zero, the voxel ID is set to 0. Otherwise, we find the traversal
    // point of the ray and the sphere center with the projected circle given by the entry_radius.
    inline int initializeAngularVoxelID(const SphericalVoxelGrid &grid, std::size_t number_of_sections,
                                        const FreeVec3 &ray_sphere, const std::vector<LineSegment> &angular_max,
                                        double ray_sphere_2, double grid_sphere_2, double entry_radius) noexcept {
        if (number_of_sections == 1) { return 0; }
        const double SED = ray_sphere.x() * ray_sphere.x() + ray_sphere_2 * ray_sphere_2;
        if (svr::isEqual(SED, 0.0)) { return 0; }
        const double r = entry_radius / std::sqrt(SED);
        const double p1 = grid.sphereCenter().x() - ray_sphere.x() * r;
        const double p2 = grid_sphere_2 - ray_sphere_2 * r;
        return calculateAngularVoxelIDFromPoints(angular_max, p1, p2);
    }

    // Determines whether a radial hit occurs for the given ray. A radial hit is considered an intersection with
    // the ray and a radial section. To determine line-sphere intersection, this follows closely the mathematics
    // presented in: http://cas.xav.free.fr/Graphics%20Gems%204%20-%20Paul%20S.%20Heckbert.pdf
    inline HitParameters radialHit(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                   RadialHitMetadata &rh_metadata, int current_radial_voxel,
                                   double v, double rsvd_minus_v_squared, double t, double t_end) noexcept {
        if (rh_metadata.radialStepTransitionHasOccurred()) {
            const double d_b = std::sqrt(grid.deltaRadiiSquared(current_radial_voxel - 1) - rsvd_minus_v_squared);
            const double intersection_t = ray.timeOfIntersectionAt(v + d_b);
            if (intersection_t < t_end) {
                return {.tMax=intersection_t, .tStep=-1, .within_bounds=true };
            }
        } else {
            const std::size_t previous_idx = std::min(static_cast<std::size_t>(current_radial_voxel),
                                                      grid.numRadialSections() - 1);
            rh_metadata.updatePreviousRadialVoxel(current_radial_voxel);
            const double r_a = grid.deltaRadiiSquared(previous_idx -
                                                      (grid.deltaRadiiSquared(previous_idx) < rsvd_minus_v_squared));
            const double d_a = std::sqrt(r_a - rsvd_minus_v_squared);
            const double t_entrance = ray.timeOfIntersectionAt(v - d_a);
            const double t_exit = ray.timeOfIntersectionAt(v + d_a);

            const double t_entrance_gt_t = t_entrance > t;
            if (t_entrance_gt_t && t_entrance == t_exit) {
                // Tangential hit.
                rh_metadata.isRadialStepTransition(true);
                return {.tMax=t_entrance, .tStep=0, .within_bounds=true};
            }
            if (t_entrance_gt_t && t_entrance < t_end) {
                return {.tMax=t_entrance, .tStep=1, .within_bounds=true};
            }
            if (t_exit < t_end) {
                // t_exit is the "further" point of intersection of the current sphere.
                // Since t_entrance is not within our time bounds, it must be true that this is a radial transition.
                rh_metadata.isRadialStepTransition(true);
                return {.tMax=t_exit, .tStep=-1, .within_bounds=true};
            }
        }
        // There does not exist an intersection time X such that t < X < t_end.
        return {.tMax=std::numeric_limits<double>::max(), .tStep=0, .within_bounds=false};
    }

    // A generalized version of the latter half of the polar and azimuthal hit parameters. Since the only difference
    // is the 2-d plane for which they exist in, this portion can be generalized to a single function.
    // The calculations presented below follow closely the works of [Foley et al, 1996], [O'Rourke, 1998].
    // Reference: http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
    HitParameters angularHit(const svr::SphericalVoxelGrid &grid, const Ray &ray, double perp_uv_min,
                             double perp_uv_max, double perp_uw_min, double perp_uw_max, double perp_vw_min,
                             double perp_vw_max, const RaySegment &ray_segment,
                             const std::array<double, 2> &collinear_times, double t, double t_end,
                             double ray_direction_2, double sphere_center_2,
                             const std::vector<svr::LineSegment> &P_max, int current_voxel) noexcept {
        const bool is_parallel_min = svr::isEqual(perp_uv_min, 0.0);
        const bool is_collinear_min = is_parallel_min && svr::isEqual(perp_uw_min, 0.0) &&
                                      svr::isEqual(perp_vw_min, 0.0);
        const bool is_parallel_max = svr::isEqual(perp_uv_max, 0.0);
        const bool is_collinear_max = is_parallel_max && svr::isEqual(perp_uw_max, 0.0) &&
                                      svr::isEqual(perp_vw_max, 0.0);
        double a, b;
        double t_min = collinear_times[is_collinear_min];
        bool is_intersect_min = false;
        if (!is_parallel_min) {
            const double inv_perp_uv_min = 1.0 / perp_uv_min;
            a = perp_vw_min * inv_perp_uv_min;
            b = perp_uw_min * inv_perp_uv_min;
            if (!((svr::lessThan(a, 0.0) || svr::lessThan(1.0, a)) || svr::lessThan(b, 0.0) || svr::lessThan(1.0, b))) {
                is_intersect_min = true;
                t_min = ray_segment.intersectionTimeAt(b, ray);
            }
        }
        double t_max = collinear_times[is_collinear_max];
        bool is_intersect_max = false;
        if (!is_parallel_max) {
            const double inv_perp_uv_max = 1.0 / perp_uv_max;
            a = perp_vw_max * inv_perp_uv_max;
            b = perp_uw_max * inv_perp_uv_max;
            if (!((svr::lessThan(a, 0.0) || svr::lessThan(1.0, a)) || svr::lessThan(b, 0.0) || svr::lessThan(1.0, b))) {
                is_intersect_max = true;
                t_max = ray_segment.intersectionTimeAt(b, ray);
            }
        }

        const bool t_max_within_bounds = svr::lessThan(t, t_max) && svr::lessThan(t_max, t_end);
        const bool t_min_within_bounds = svr::lessThan(t, t_min) && svr::lessThan(t_min, t_end);
        if (!t_max_within_bounds && !t_min_within_bounds) {
            return {.tMax = std::numeric_limits<double>::max(), .tStep = 0, .within_bounds = false};
        }

        if (is_intersect_max && !is_intersect_min && !is_collinear_min) {
            return {.tMax = t_max, .tStep = 1, .within_bounds = t_max_within_bounds};
        }
        if (is_intersect_min && !is_intersect_max && !is_collinear_max) {
            return {.tMax = t_min, .tStep = -1, .within_bounds = t_min_within_bounds};
        }
        if ((is_intersect_min && is_intersect_max) ||
            (is_intersect_min && is_collinear_max) ||
            (is_intersect_max && is_collinear_min)) {
            if (svr::isEqual(t_min, t_max)) {
                const double perturbed_t = 0.1;
                a = -ray.direction().x() * perturbed_t;
                b = -ray_direction_2 * perturbed_t;
                const double max_radius_over_plane_length = grid.sphereMaxRadius() / std::sqrt(a * a + b * b);
                const double p1 = grid.sphereCenter().x() - max_radius_over_plane_length * a;
                const double p2 = sphere_center_2 - max_radius_over_plane_length * b;
                const int next_step = std::abs(current_voxel - calculateAngularVoxelIDFromPoints(P_max, p1, p2));
                return {.tMax = t_max,
                        .tStep = (svr::lessThan(ray_direction_2, 0.0) || svr::lessThan(ray.direction().x(), 0.0)) ?
                                 next_step : -next_step,
                        .within_bounds = t_min_within_bounds
                };
            }
            if (t_min_within_bounds && (svr::lessThan(t_min, t_max) || svr::isEqual(t, t_max))) {
                return {.tMax = t_min, .tStep = -1, .within_bounds = true};
            }
            if (t_max_within_bounds && (svr::lessThan(t_max, t_min) || svr::isEqual(t, t_min))) {
                return {.tMax = t_max, .tStep = 1, .within_bounds = true};
            }
        }
        return {.tMax = std::numeric_limits<double>::max(), .tStep = 0, .within_bounds = false};
    }

    // Determines whether a polar hit occurs for the given ray. A polar hit is considered an intersection with
    // the ray and a polar section. The polar sections live in the XY plane.
    inline HitParameters polarHit(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                  const RaySegment &ray_segment, const std::array<double, 2> &collinear_times,
                                  int current_polar_voxel, double t, double t_end) noexcept {
        // Calculate the voxel boundary vectors.
        const FreeVec3 p_one(grid.pMaxPolar(current_polar_voxel).P1,
                             grid.pMaxPolar(current_polar_voxel).P2, 0.0);
        const FreeVec3 p_two(grid.pMaxPolar(current_polar_voxel + 1).P1,
                             grid.pMaxPolar(current_polar_voxel + 1).P2, 0.0);
        const BoundVec3 *u_min = &grid.centerToPolarBound(current_polar_voxel);
        const BoundVec3 *u_max = &grid.centerToPolarBound(current_polar_voxel + 1);
        const FreeVec3 w_min = p_one - FreeVec3(ray_segment.P1());
        const FreeVec3 w_max = p_two - FreeVec3(ray_segment.P1());
        const double perp_uv_min = u_min->x() * ray_segment.vector().y() - u_min->y() * ray_segment.vector().x();
        const double perp_uv_max = u_max->x() * ray_segment.vector().y() - u_max->y() * ray_segment.vector().x();
        const double perp_uw_min = u_min->x() * w_min.y() - u_min->y() * w_min.x();
        const double perp_uw_max = u_max->x() * w_max.y() - u_max->y() * w_max.x();
        const double perp_vw_min = ray_segment.vector().x() * w_min.y() - ray_segment.vector().y() * w_min.x();
        const double perp_vw_max = ray_segment.vector().x() * w_max.y() - ray_segment.vector().y() * w_max.x();
        return angularHit(grid, ray, perp_uv_min, perp_uv_max, perp_uw_min,
                          perp_uw_max, perp_vw_min, perp_vw_max, ray_segment,
                          collinear_times, t, t_end, ray.direction().y(),
                          grid.sphereCenter().y(), grid.pMaxPolar(),
                          current_polar_voxel);
    }

    // Determines whether an azimuthal hit occurs for the given ray. An azimuthal hit is
    // considered an intersection with the ray and an azimuthal section.
    // The azimuthal sections live in the XZ plane.
    inline HitParameters azimuthalHit(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                      const RaySegment &ray_segment, const std::array<double, 2> &collinear_times,
                                      int current_azimuthal_voxel, double t, double t_end) noexcept {
        // Calculate the voxel boundary vectors.
        const FreeVec3 p_one(grid.pMaxAzimuthal(current_azimuthal_voxel).P1, 0.0,
                             grid.pMaxAzimuthal(current_azimuthal_voxel).P2);
        const FreeVec3 p_two(grid.pMaxAzimuthal(current_azimuthal_voxel + 1).P1, 0.0,
                             grid.pMaxAzimuthal(current_azimuthal_voxel + 1).P2);
        const BoundVec3 *u_min = &grid.centerToAzimuthalBound(current_azimuthal_voxel);
        const BoundVec3 *u_max = &grid.centerToAzimuthalBound(current_azimuthal_voxel + 1);
        const FreeVec3 w_min = p_one - FreeVec3(ray_segment.P1());
        const FreeVec3 w_max = p_two - FreeVec3(ray_segment.P1());
        const double perp_uv_min = u_min->x() * ray_segment.vector().z() - u_min->z() * ray_segment.vector().x();
        const double perp_uv_max = u_max->x() * ray_segment.vector().z() - u_max->z() * ray_segment.vector().x();
        const double perp_uw_min = u_min->x() * w_min.z() - u_min->z() * w_min.x();
        const double perp_uw_max = u_max->x() * w_max.z() - u_max->z() * w_max.x();
        const double perp_vw_min = ray_segment.vector().x() * w_min.z() - ray_segment.vector().z() * w_min.x();
        const double perp_vw_max = ray_segment.vector().x() * w_max.z() - ray_segment.vector().z() * w_max.x();
        return angularHit(grid, ray, perp_uv_min, perp_uv_max, perp_uw_min,
                          perp_uw_max, perp_vw_min, perp_vw_max, ray_segment,
                          collinear_times, t, t_end, ray.direction().z(),
                          grid.sphereCenter().z(), grid.pMaxAzimuthal(),
                          current_azimuthal_voxel);
    }

    // Calculates the voxel(s) with the minimal tMax for the next intersection. Since t is being updated
    // with each interval of the algorithm, this must check the following cases:
    // 1. tMaxR is the minimum.
    // 2. tMaxTheta is the minimum.
    // 3. tMaxPhi is the minimum.
    // 4. tMaxR, tMaxTheta, tMaxPhi equal intersection.
    // 5. tMaxR, tMaxTheta equal intersection.
    // 6. tMaxR, tMaxPhi equal intersection.
    // 7. tMaxTheta, tMaxPhi equal intersection.
    //
    // For each case, the following must hold: t < tMax < t_end
    inline VoxelIntersectionType minimumIntersection(const HitParameters &rad_params,
                                                     const HitParameters &ang_params,
                                                     const HitParameters &azi_params) noexcept {
        if (rad_params.within_bounds && svr::lessThan(rad_params.tMax, ang_params.tMax)
            && svr::lessThan(rad_params.tMax, azi_params.tMax)) {
            return VoxelIntersectionType::Radial;
        }
        if (ang_params.within_bounds && svr::lessThan(ang_params.tMax, rad_params.tMax)
            && svr::lessThan(ang_params.tMax, azi_params.tMax)) {
            return VoxelIntersectionType::Polar;
        }
        if (azi_params.within_bounds && svr::lessThan(azi_params.tMax, ang_params.tMax)
            && svr::lessThan(azi_params.tMax, rad_params.tMax)) {
            return VoxelIntersectionType::Azimuthal;
        }
        if (rad_params.within_bounds && svr::isEqual(rad_params.tMax, ang_params.tMax)
            && isEqual(rad_params.tMax, azi_params.tMax)) {
            return VoxelIntersectionType::RadialPolarAzimuthal;
        }
        if (azi_params.within_bounds && svr::isEqual(azi_params.tMax, ang_params.tMax)) {
            return VoxelIntersectionType::PolarAzimuthal;
        }
        if (rad_params.within_bounds && svr::isEqual(ang_params.tMax, rad_params.tMax)) {
            return VoxelIntersectionType::RadialPolar;
        }
        if (rad_params.within_bounds && svr::isEqual(rad_params.tMax, azi_params.tMax)) {
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
        if (grid.numPolarSections() == grid.numAzimuthalSections()) {
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
                                                         double t_begin, double t_end) noexcept {
        const FreeVec3 rsv_begin = grid.sphereCenter() - ray.pointAtParameter(t_begin); // Ray Sphere Vector.
        const double SED_from_center = rsv_begin.squared_length();
        std::size_t radial_entrance_voxel = 0;
        const double max_radius_squared = grid.deltaRadiiSquared(0);

        while (SED_from_center < grid.deltaRadiiSquared(radial_entrance_voxel)) { ++radial_entrance_voxel; }
        const bool ray_origin_is_outside_grid = (radial_entrance_voxel == 0);

        const std::size_t vector_index = radial_entrance_voxel - !ray_origin_is_outside_grid;
        const double entry_radius = grid.deltaRadii(vector_index);
        const double entry_radius_squared = grid.deltaRadiiSquared(vector_index);

        const FreeVec3 rsv = t_begin == 0.0 ? rsv_begin : grid.sphereCenter() - ray.pointAtParameter(0.0);
        const double rsvd = rsv.dot(rsv); // Ray Sphere Vector Dot product at time zero.
        const double v = rsv.dot(ray.unitDirection().to_free());
        const double rsvd_minus_v_squared = rsvd - v * v;

        if (entry_radius_squared <= rsvd_minus_v_squared) { return {}; }
        const double d = std::sqrt(entry_radius_squared - rsvd_minus_v_squared);

        const double t_sphere_entrance = ray.timeOfIntersectionAt(v - d);
        const double t_sphere_exit = ray.timeOfIntersectionAt(v + d);

        if (t_sphere_entrance < t_begin && t_sphere_exit < t_begin) { return {}; }
        int current_radial_voxel = radial_entrance_voxel + ray_origin_is_outside_grid;

        std::vector<svr::LineSegment> P_polar(grid.numPolarSections() + 1);
        std::vector<svr::LineSegment> P_azimuthal(grid.numAzimuthalSections() + 1);
        initializeVoxelBoundarySegments(P_polar, P_azimuthal, ray_origin_is_outside_grid, grid, entry_radius);

        const FreeVec3 ray_sphere = ray_origin_is_outside_grid ?
                                    grid.sphereCenter() - ray.pointAtParameter(t_sphere_entrance) :
                                    SED_from_center == 0.0 ?
                                    grid.sphereCenter() - ray.pointAtParameter(t_begin + 0.1) : rsv_begin;

        int current_polar_voxel = initializeAngularVoxelID(grid, grid.numPolarSections(), ray_sphere, P_polar,
                                                           ray_sphere.y(), grid.sphereCenter().y(), entry_radius);
        if (static_cast<std::size_t>(current_polar_voxel) == grid.numPolarSections()) { return {}; }

        int current_azimuthal_voxel = initializeAngularVoxelID(grid, grid.numAzimuthalSections(), ray_sphere,
                                                               P_azimuthal, ray_sphere.z(), grid.sphereCenter().z(),
                                                               entry_radius);
        if (static_cast<std::size_t>(current_azimuthal_voxel) == grid.numAzimuthalSections()) { return {}; }

        std::vector<svr::SphericalVoxel> voxels;
        voxels.reserve(grid.numRadialSections() + grid.numPolarSections() + grid.numAzimuthalSections());
        voxels.push_back({.radial=current_radial_voxel,
                          .polar=current_polar_voxel,
                          .azimuthal=current_azimuthal_voxel});
        double t;
        if (ray_origin_is_outside_grid) {
            t = t_sphere_entrance;
            t_end = std::min(t_end, t_sphere_exit);
        } else {
            t = t_begin;
            const double max_d = std::sqrt(max_radius_squared - rsvd_minus_v_squared);
            t_end = std::min(t_end, ray.timeOfIntersectionAt(v + max_d));
        }

        // Initialize the time in case of collinear min or collinear max for angular plane hits.
        // In the case where the hit is not collinear, a time of 0.0 is inputted.
        const std::array<double, 2> collinear_times = {0.0, ray.timeOfIntersectionAt(grid.sphereCenter())};

        RadialHitMetadata rh_metadata;
        rh_metadata.updatePreviousRadialVoxel(current_radial_voxel);
        RaySegment ray_segment(t_end, ray);

        while (true) {
            const auto radial = radialHit(ray, grid, rh_metadata, current_radial_voxel,
                                          v, rsvd_minus_v_squared, t, t_end);
            if (current_radial_voxel + radial.tStep == 0) { return voxels; }
            ray_segment.updateAtTime(t, ray);
            const auto polar = polarHit(ray, grid, ray_segment, collinear_times,
                                        current_polar_voxel, t, t_end);
            const auto azimuthal = azimuthalHit(ray, grid, ray_segment, collinear_times,
                                                current_azimuthal_voxel, t, t_end);
            const auto voxel_intersection = minimumIntersection(radial, polar, azimuthal);
            switch (voxel_intersection) {
                case Radial: {
                    t = radial.tMax;
                    current_radial_voxel += radial.tStep;
                    if (rh_metadata.previousRadialVoxel() == current_radial_voxel) { continue; }
                    break;
                }
                case Polar: {
                    t = polar.tMax;
                    current_polar_voxel = (current_polar_voxel + polar.tStep) % grid.numPolarSections();
                    break;
                }
                case Azimuthal: {
                    t = azimuthal.tMax;
                    current_azimuthal_voxel = (current_azimuthal_voxel + azimuthal.tStep) % grid.numAzimuthalSections();
                    break;
                }
                case RadialPolar: {
                    t = radial.tMax;
                    current_radial_voxel += radial.tStep;
                    current_polar_voxel = (current_polar_voxel + polar.tStep) % grid.numPolarSections();
                    break;
                }
                case RadialAzimuthal: {
                    t = radial.tMax;
                    current_radial_voxel += radial.tStep;
                    current_azimuthal_voxel = (current_azimuthal_voxel + azimuthal.tStep) % grid.numAzimuthalSections();
                    break;
                }
                case PolarAzimuthal: {
                    t = polar.tMax;
                    current_polar_voxel = (current_polar_voxel + polar.tStep) % grid.numPolarSections();
                    current_azimuthal_voxel = (current_azimuthal_voxel + azimuthal.tStep) % grid.numAzimuthalSections();
                    break;
                }
                case RadialPolarAzimuthal: {
                    t = radial.tMax;
                    current_radial_voxel += radial.tStep;
                    current_polar_voxel = (current_polar_voxel + polar.tStep) % grid.numPolarSections();
                    current_azimuthal_voxel = (current_azimuthal_voxel + azimuthal.tStep) % grid.numAzimuthalSections();
                    break;
                }
                case None: { return voxels; }
            }
            voxels.push_back({.radial=current_radial_voxel,
                              .polar=current_polar_voxel,
                              .azimuthal=current_azimuthal_voxel});
        }
    }

    std::vector<svr::SphericalVoxel> walkSphericalVolume(double *ray_origin, double *ray_direction,
                                                         double *min_bound, double *max_bound,
                                                         std::size_t num_radial_voxels,
                                                         std::size_t num_polar_voxels,
                                                         std::size_t num_azimuthal_voxels,
                                                         double *sphere_center, double t_begin, double t_end) noexcept {
        return svr::walkSphericalVolume(Ray(BoundVec3(ray_origin[0], ray_origin[1], ray_origin[2]),
                                            FreeVec3(ray_direction[0], ray_direction[1], ray_direction[2])),
                                        svr::SphericalVoxelGrid(svr::SphereBound{.radial=min_bound[0],
                                                                        .polar=min_bound[1],
                                                                        .azimuthal=min_bound[2]},
                                                                svr::SphereBound{.radial=max_bound[0],
                                                                        .polar=max_bound[1],
                                                                        .azimuthal=max_bound[2]},
                                                                num_radial_voxels,
                                                                num_polar_voxels,
                                                                num_azimuthal_voxels,
                                                                BoundVec3(sphere_center[0],
                                                                          sphere_center[1],
                                                                          sphere_center[2])), t_begin, t_end);
    }

} // namespace svr