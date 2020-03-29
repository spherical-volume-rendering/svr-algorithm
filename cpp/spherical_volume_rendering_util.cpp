#include "spherical_volume_rendering_util.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

// Takes the absolute value of 'value' and determines if it is less than the tolerance.
// In this case, we assume tolerance is an arbitrary small value.
// This is used when floating point mathematics carries rounding error.
inline bool isNearZero(double value, double tolerance) {
    return std::abs(value) < tolerance;
}

// Produces the corresponding angular or azimuthal voxel ID given two directions.
// This is done using the four-quadrant inverse tangent of a two-dimensional plane. From this,
// it will return signed values in the counterclockwise direction of the angles in the Euclidean plane between
// the ray origin on the closed interval [-pi, pi].
//
// Note: With the current variable naming (x,y), we assume the point is on the XY plane.
//       This function will work with any 2-dimensional plane.
//       If you want a point on the XZ plane, then y is the z-coordinate, x is the x-coordinate.
inline size_t getInitialVoxelID(double y, double x, size_t num_sections) {
    const size_t current_voxel_ID = std::floor(std::atan2(y, x) * num_sections / ( 2 * M_PI));
    return (current_voxel_ID + num_sections) % num_sections;
}

// The parameters returned by radialHit().
struct RadialHitParameters {
    // The time at which a hit occurs for the ray at the next point of intersection with a radial section.
    // This is always calculated from the ray origin.
    double tMaxR;
    // The voxel traversal value of a radial step: 0, +1, -1. This is added to the current radial voxel.
    size_t tStepR;
    // Determine whether the current voxel traversal was a 'transition'. This is necessary to determine when
    // The radial steps should go from negative to positive or vice-versa, since the radial voxels go from 1..N, where N
    // is the number of radial sections.
    bool previous_transition_flag;
};

// The parameters returned by angularHit().
struct AngularHitParameters {
    // The time at which a hit occurs for the ray at the next point of intersection with an angular section.
    // This is always calculated from the ray origin.
    double tMaxTheta;
    // The voxel traversal value of an angular step: 0, +1, -1. This is added to the current angular voxel.
    size_t tStepTheta;
};

// The parameters returned by azimuthalHit().
struct AzimuthalHitParameters {
    // The time at which a hit occurs for the ray at the next point of intersection with an azimuthal section.
    // This is always calculated from the ray origin.
    double tMaxPhi;
    // The voxel traversal value of an azimuthal step: 0, +1, -1. This is added to the current azimuthal voxel.
    size_t tStepPhi;
};

// Determines whether a radial hit occurs for the given ray. A radial hit is considered an intersection with
// the ray and a radial section. This follows closely the mathematics presented in:
// http://cas.xav.free.fr/Graphics%20Gems%204%20-%20Paul%20S.%20Heckbert.pdf
// Input:
//    ray: The given ray to check for intersection.
//    grid: The grid that the ray is intersecting with.
//    current_voxel_ID_r: The current radial voxel ID.
//    ray_sphere_vector_dot: The dot product of the vector difference between the ray origin and the sphere origin.
//    t: The current time.
//    v: The dot product between the ray's direction and the ray_sphere_vector.
//    tol: The allowed tolerance for float point error.
//    prev_transition_flag: Determines whether the previous radial traversal was a 'transition'. A transition
//                          is defined as the change in sign of tStepR. Another way this can be determined is
//                          sequentil hits with equal radii.
//
// Returns: The corresponding radial hit parameters.
RadialHitParameters radialHit(const Ray& ray, const SphericalVoxelGrid& grid, size_t current_voxel_ID_r,
                              double ray_sphere_vector_dot, double t, double v, double tol,
                              bool previous_transition_flag) noexcept {
    const double current_radius = grid.sphereMaxRadius() - grid.deltaRadius() * (current_voxel_ID_r - 1);
    double r_a = std::max(current_radius - grid.deltaRadius(), grid.deltaRadius());

    double r_b;
    if (!previous_transition_flag) {
        // To find the next radius, we need to check the previous_transition_flag:
        // In the case that the ray has sequential hits with equal radii, e.g.
        // the innermost radial disc, this ensures that the proper radii are being checked.
        r_b = std::min(current_radius + grid.deltaRadius(), grid.sphereMaxRadius());
    } else { r_b = std::min(current_radius, grid.sphereMaxRadius());  }

    // Find the intersection times for the ray and the previous and next radial discs.
    const double ray_sphere_dot_minus_v_squared = ray_sphere_vector_dot - v * v;
    double discriminant_a = r_a * r_a - ray_sphere_dot_minus_v_squared;
    if (discriminant_a < 0) {
        r_a += grid.deltaRadius();
        discriminant_a = r_a * r_a - ray_sphere_dot_minus_v_squared;
    }
    double t1;
    double t2;
    const double d_a = std::sqrt(discriminant_a);
    t1 = ray.timeOfIntersectionAt(v - d_a);
    t2 = ray.timeOfIntersectionAt(v + d_a);
    std::vector<double> intersection_times(4);
    intersection_times[0] = t1;
    intersection_times[1] = t2;

    const double discriminant_b = r_b * r_b - ray_sphere_dot_minus_v_squared;
    if (discriminant_b >= 0) {
        const double d_b = std::sqrt(discriminant_b);
        const BoundVec3 pa = ray.origin() + ray.direction() * (v - d_b);
        const BoundVec3 pb = ray.origin() + ray.direction() * (v + d_b);
        t1 = ray.timeOfIntersectionAt(v - d_b);
        t2 = ray.timeOfIntersectionAt(v + d_b);
        intersection_times[2] = t1;
        intersection_times[3] = t2;
    }
    intersection_times.erase(std::remove_if(intersection_times.begin(), intersection_times.end(),
                                            [t](int i) { return i <= t; }), intersection_times.end());

    RadialHitParameters radial_params;
    if (intersection_times.size() >= 2 && isNearZero(intersection_times[0] - intersection_times[1], tol)) {
        // Ray is tangent to the circle, i.e. two intersection times are equal.
        radial_params.tMaxR = intersection_times[0];
        const BoundVec3 p = ray.origin() + ray.direction() * radial_params.tMaxR;
        const double r_new = std::sqrt(p.x() - grid.sphereCenter().x() + p.y() - grid.sphereCenter().y() +
                                       p.z() - grid.sphereCenter().z());
        radial_params.tStepR = 0;
        if (isNearZero(current_radius - r_new, tol)) {
            radial_params.previous_transition_flag = true;
        } else {
            radial_params.previous_transition_flag = false;
        }
    } else if (intersection_times.empty()) {
        // No intersection.
        radial_params.tMaxR = std::numeric_limits<double>::infinity();
        radial_params.tStepR = 0;
        radial_params.previous_transition_flag = false;
    } else {
        // Radial intersection.
        radial_params.tMaxR = intersection_times[0];
        const BoundVec3 p = ray.origin() + ray.direction() * radial_params.tMaxR;
        const double r_new = std::sqrt(p.x() - grid.sphereCenter().x() + p.y() - grid.sphereCenter().y() +
                                       p.z() - grid.sphereCenter().z());

        const double radial_difference = r_new - current_radius;
        const bool radial_difference_is_near_zero = isNearZero(radial_difference, tol);
        if (radial_difference_is_near_zero) {
            radial_params.previous_transition_flag = true;
        } else {
            radial_params.previous_transition_flag = false;
        }

        if (radial_difference < 0 && std::abs(radial_difference) > tol && !radial_difference_is_near_zero) {
            radial_params.tStepR = 1;
        } else {
            radial_params.tStepR = -1;
        }
    }
    return radial_params;
}

// Determines whether an angular hit occurs for the given ray. An angular hit is considered an intersection with
// the ray and an angular section. The angular sections live in the XY plane.
// Returns: the corresponding angular hit parameters.
AngularHitParameters angularHit(const Ray& ray, const SphericalVoxelGrid& grid,
        size_t current_voxel_ID_theta, double t) noexcept {
    assert(false);
}

// Determines whether an azimuthal hit occurs for the given ray. An azimuthal hit is considered an intersection with
// the ray and an azimuthal section. The azimuthal sections live in the XZ plane.
// Returns: the corresponding azimuthal hit parameters.
AzimuthalHitParameters azimuthalHit(const Ray& ray, const SphericalVoxelGrid& grid,
        size_t current_phi_ID_theta, double t) noexcept {
    assert(false);
}

std::vector<SphericalVoxel>
sphericalCoordinateVoxelTraversal(const Ray &ray, const SphericalVoxelGrid &grid, double t_begin,
                                       double t_end, double tol) noexcept {
    std::vector<SphericalVoxel> voxels;
    voxels.reserve(grid.numRadialVoxels() + grid.numAngularVoxels() + grid.numAzimuthalVoxels());

    /* INITIALIZATION PHASE */

    // Determine ray location at t_begin.
    const BoundVec3 point_at_t_begin = ray.pointAtParameter(t_begin);
    const FreeVec3 ray_sphere_vector = grid.sphereCenter() - point_at_t_begin;
    double current_r = grid.deltaRadius();

    const double ray_sphere_vector_dot = ray_sphere_vector.dot(ray_sphere_vector);
    while (ray_sphere_vector_dot > (current_r * current_r) && current_r < grid.sphereMaxRadius()) {
        current_r += grid.deltaRadius();
    }

    // Find the intersection times for the ray and the radial shell containing the parameter point at t_begin.
    // This will determine if the ray intersects the sphere.
    const double v = ray_sphere_vector.dot(ray.direction().to_free());
    const double discriminant = (current_r * current_r) - ray_sphere_vector_dot - (v * v);
    if (discriminant <= 0) { return voxels; }
    const double d = std::sqrt(discriminant);

    // Calculate the time of entrance and exit of the ray.
    // Need to use a non-zero direction to determine this.
    double t1 = ray.timeOfIntersectionAt(v - d);
    double t2 = ray.timeOfIntersectionAt(v + d);

    if ((t1 < t_begin && t2 < t_begin) || isNearZero(t1 - t2, tol)) {
        // Case 1: No intersection.
        // Case 2: Tangent hit.
        return voxels;
    }
    size_t current_voxel_ID_r = 1 + (grid.sphereMaxRadius() - current_r) / grid.deltaRadius();

    size_t current_voxel_ID_theta;
    size_t current_voxel_ID_phi;
    if (isNearZero(ray.origin().x() - grid.sphereCenter().x(), tol) &&
        isNearZero(ray.origin().y() - grid.sphereCenter().y(), tol) &&
        isNearZero(ray.origin().z() - grid.sphereCenter().z(), tol)) {
        // If the ray starts at the sphere's center, we need to perturb slightly along
        // the path to determine the correct angular and azimuthal voxel.
        const double perturbed_t = 0.1;
        const double perturbed_x = point_at_t_begin.x() + ray.direction().x() * perturbed_t;
        const double perturbed_y = point_at_t_begin.y() + ray.direction().y() * perturbed_t;
        const double perturbed_z = point_at_t_begin.z() + ray.direction().z() * perturbed_t;
        current_voxel_ID_theta = getInitialVoxelID(perturbed_y, perturbed_x, grid.numAngularVoxels());
        current_voxel_ID_phi = getInitialVoxelID(perturbed_z, perturbed_x, grid.numAzimuthalVoxels());
    } else {
        current_voxel_ID_theta = getInitialVoxelID(point_at_t_begin.y() - grid.sphereCenter().y(),
                                                       point_at_t_begin.x() - grid.sphereCenter().x(),
                                                       grid.numAngularVoxels());
        current_voxel_ID_phi = getInitialVoxelID(point_at_t_begin.z() - grid.sphereCenter().z(),
                                                   point_at_t_begin.x() - grid.sphereCenter().x(),
                                                   grid.numAzimuthalVoxels());
    }

    voxels.push_back({.radial_voxel=current_voxel_ID_r,
                      .angular_voxel=current_voxel_ID_theta,
                      .azimuthal_voxel=current_voxel_ID_phi});

    // Find the maximum time the ray will be in the grid.
    const double max_discriminant = grid.sphereMaxRadius() * grid.sphereMaxRadius() - (ray_sphere_vector_dot - v * v);
    const double max_d = std::sqrt(max_discriminant);
    const BoundVec3 pa_max = ray.origin() + ray.direction() * (v - max_d);
    const BoundVec3 pb_max = ray.origin() + ray.direction() * (v + max_d);
    t1 = ray.timeOfIntersectionAt(v - max_d);
    t2 = ray.timeOfIntersectionAt(v + max_d);
    const double t_grid_end = std::min(t1, t2);

    /* TRAVERSAL PHASE */
    double t = t_begin;
    t_end = std::min(t_grid_end, t_end);
    bool previous_transition_flag = false;

    while (t < t_end) {
        const RadialHitParameters radial_params = radialHit(ray, grid, current_voxel_ID_r,
                ray_sphere_vector_dot, t, v, tol, previous_transition_flag);
        previous_transition_flag = radial_params.previous_transition_flag;
        const AngularHitParameters angular_params = angularHit(ray, grid, current_voxel_ID_theta, t);
        const AzimuthalHitParameters azimuthal_params = azimuthalHit(ray, grid, current_voxel_ID_phi, t);

        if (radial_params.tMaxR <= angular_params.tMaxTheta &&
            radial_params.tMaxR <= azimuthal_params.tMaxPhi &&
            t < radial_params.tMaxR && radial_params.tMaxR < t_end &&
            current_voxel_ID_r + radial_params.tStepR != 0) {
            // tMaxR is the minimum, the next radial step is within bounds (t, t_end),
            // and the next step is not a radial exit.
            t = radial_params.tMaxR;
            current_voxel_ID_r += radial_params.tStepR;
        } else if (angular_params.tMaxTheta <= azimuthal_params.tMaxPhi &&
                  t < angular_params.tMaxTheta && angular_params.tMaxTheta <= t_end) {
            // tMaxTheta is the minimum and the next angular step is within bounds (t, t_end).
            t = angular_params.tMaxTheta;
            current_voxel_ID_theta += angular_params.tStepTheta % grid.numAngularVoxels();
        } else if (t < azimuthal_params.tMaxPhi && azimuthal_params.tMaxPhi < t_end) {
            // tMaxPhi is the minimum and the next azimuthal step is within bounds (t, t_end).
            t = azimuthal_params.tMaxPhi;
            current_voxel_ID_phi += azimuthal_params.tStepPhi % grid.numAzimuthalVoxels();
        } else {
            // No hits are within the bounds (t, t_end).
            return voxels;
        }
        voxels.push_back({.radial_voxel=current_voxel_ID_r,
                          .angular_voxel=current_voxel_ID_theta,
                          .azimuthal_voxel=current_voxel_ID_phi});
    }
    return voxels;
}

std::vector<SphericalVoxel> sphericalCoordinateVoxelTraversalCy(double* ray_origin, double* ray_direction,
                                                                double* min_bound, double* max_bound,
                                                                size_t num_radial_voxels, size_t num_angular_voxels,
                                                                size_t num_azimuthal_voxels, double* sphere_center,
                                                                double sphere_max_radius, double t_begin, double t_end,
                                                                double tol) noexcept {
    const BoundVec3 ray_origin_t (ray_origin[0], ray_origin[1], ray_origin[2]);
    const UnitVec3 ray_direction_t (ray_direction[0], ray_direction[1], ray_direction[2]);
    const Ray ray(ray_origin_t, ray_direction_t);

    const BoundVec3 max_bound_t(max_bound[0], max_bound[1], max_bound[2]);
    const BoundVec3 min_bound_t(min_bound[0], min_bound[1], min_bound[2]);
    const BoundVec3 sphere_center_t(sphere_center[0], sphere_center[1], sphere_center[2]);
    const SphericalVoxelGrid grid(min_bound_t, max_bound_t, num_radial_voxels, num_angular_voxels,
                                  num_azimuthal_voxels, sphere_center_t, sphere_max_radius);

    return sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end, tol);
}