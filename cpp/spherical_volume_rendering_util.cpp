#include "spherical_volume_rendering_util.h"
#include <vector>
#include <cmath>

// Takes the absolute value of 'value' and determines if it is less than the tolerance.
// In this case, we assume tolerance is an arbitrary small value.
// This is used when floating point mathematics carries rounding error.
inline bool isNearZero(double value, double tolerance) {
    return std::abs(value) < tolerance;
}

// Produces the corresponding angular or azimuthal voxel ID given two directions.
// TODO(cgyurgyik): Update documentation, provide link to better description.
inline size_t getInitialVoxelID(double y, double x, size_t num_sections) {
    const size_t current_voxel_ID = std::floor(std::atan2(y, x) * num_sections / ( 2 * M_PI));
    return (current_voxel_ID + num_sections) % num_sections;
}

struct RadialHitParameters {
    double tMaxR;
    size_t tStepR;
    bool previous_transition_flag;
};

struct AngularHitParameters {
    double tMaxTheta;
    size_t tStepTheta;
};

struct AzimuthalHitParameters {
    double tMaxPhi;
    size_t tStepPhi;
};

// TODO(cgyurgyik): Documentation, implementation.
RadialHitParameters radialHit(const Ray& ray, const SphericalVoxelGrid& grid, size_t current_voxel_ID_r,
                              const FreeVec3& ray_sphere_vector, double t, double v,
                              bool previous_transition_flag) noexcept {
    assert(false);
}

// TODO(cgyurgyik): Documentation, implementation.
AngularHitParameters angularHit(const Ray& ray, const SphericalVoxelGrid& grid,
        size_t current_voxel_ID_theta, double t) noexcept {
    assert(false);
}

// TODO(cgyurgyik): Documentation, implementation.
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
    const BoundVec3 point_at_t_begin = ray.point_at_parameter(t_begin);
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
    const BoundVec3 pa = ray.origin() + ray.direction() * (v - d);
    const BoundVec3 pb = ray.origin() + ray.direction() * (v + d);

    // Calculate the time of entrance and exit of the ray.
    // Need to use a non-zero direction to determine this.
    double t1;
    double t2;
    if (std::abs(ray.direction().y()) > tol) {
        t1 = (pa.y() - ray.origin().y()) * ray.inv_direction().y();
        t2 = (pb.y() - ray.origin().y()) * ray.inv_direction().y();
    } else if (std::abs(ray.direction().x()) > tol) {
        t1 = (pa.x() - ray.origin().x()) * ray.inv_direction().x();
        t2 = (pb.x() - ray.origin().x()) * ray.inv_direction().x();
    } else {
        t1 = (pa.z() - ray.origin().z()) * ray.inv_direction().z();
        t2 = (pb.z() - ray.origin().z()) * ray.inv_direction().z();
    }

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
    if (std::abs(ray.direction().y()) > tol) {
        t1 = (pa_max.y() - ray.origin().y()) * ray.inv_direction().y();
        t2 = (pb_max.y() - ray.origin().y()) * ray.inv_direction().y();
    } else if (std::abs(ray.direction().x()) > tol) {
        t1 = (pa_max.x() - ray.origin().x()) * ray.inv_direction().x();
        t2 = (pb_max.x() - ray.origin().x()) * ray.inv_direction().x();
    } else {
        t1 = (pa_max.z() - ray.origin().z()) * ray.inv_direction().z();
        t2 = (pb_max.z() - ray.origin().z()) * ray.inv_direction().z();
    }
    const double t_grid_end = std::min(t1, t2);

    /* TRAVERSAL PHASE */
    double t = t_begin;
    t_end = std::min(t_grid_end, t_end);
    bool previous_transition_flag = false;

    while (t < t_end) {
        const RadialHitParameters radial_params = radialHit(ray, grid, current_voxel_ID_r,
                ray_sphere_vector, t, v, previous_transition_flag);
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

    std::vector<SphericalVoxel> voxels;
    voxels.reserve(grid.numRadialVoxels() + grid.numAngularVoxels() + grid.numAzimuthalVoxels());

    /* INITIALIZATION PHASE */

    // Determine ray location at t_begin.
    const BoundVec3 point_at_t_begin = ray.point_at_parameter(t_begin);
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
    const BoundVec3 pa = ray.origin() + ray.direction() * (v - d);
    const BoundVec3 pb = ray.origin() + ray.direction() * (v + d);

    // Calculate the time of entrance and exit of the ray.
    // Need to use a non-zero direction to determine this.
    double t1;
    double t2;
    if (std::abs(ray.direction().y()) > tol) {
        t1 = (pa.y() - ray.origin().y()) * ray.inv_direction().y();
        t2 = (pb.y() - ray.origin().y()) * ray.inv_direction().y();
    } else if (std::abs(ray.direction().x()) > tol) {
        t1 = (pa.x() - ray.origin().x()) * ray.inv_direction().x();
        t2 = (pb.x() - ray.origin().x()) * ray.inv_direction().x();
    } else {
        t1 = (pa.z() - ray.origin().z()) * ray.inv_direction().z();
        t2 = (pb.z() - ray.origin().z()) * ray.inv_direction().z();
    }

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
    if (std::abs(ray.direction().y()) > tol) {
        t1 = (pa_max.y() - ray.origin().y()) * ray.inv_direction().y();
        t2 = (pb_max.y() - ray.origin().y()) * ray.inv_direction().y();
    } else if (std::abs(ray.direction().x()) > tol) {
        t1 = (pa_max.x() - ray.origin().x()) * ray.inv_direction().x();
        t2 = (pb_max.x() - ray.origin().x()) * ray.inv_direction().x();
    } else {
        t1 = (pa_max.z() - ray.origin().z()) * ray.inv_direction().z();
        t2 = (pb_max.z() - ray.origin().z()) * ray.inv_direction().z();
    }
    const double t_grid_end = std::min(t1, t2);

    /* TRAVERSAL PHASE */
    double t = t_begin;
    t_end = std::min(t_grid_end, t_end);
    bool previous_transition_flag = false;

    while (t < t_end) {
        const RadialHitParameters radial_params = radialHit(ray, grid, current_voxel_ID_r,
                                                            ray_sphere_vector, t, v, previous_transition_flag);
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