#include "spherical_volume_rendering_util.h"
#include <vector>
#include <cmath>

// Takes the absolute value of 'value' and determines if it is less than the tolerance.
// In this case, we assume tolerance is an arbitrary small value.
// This is used when floating point mathematics carries rounding error.
[[nodiscard]] inline bool isNearZero(double value, double tolerance) {
    return std::abs(value) < tolerance;
}

// Produces the corresponding angular or azimuthal voxel ID given two directions.
// TODO(cgyurgyik): Update documentation, provide link to better description.
[[nodiscard]] inline size_t getInitialVoxelID(double y, double x, size_t num_sections) {
    const size_t current_voxel_ID = std::floor(std::atan2(y, x) * num_sections / ( 2 * M_PI));
    if (current_voxel_ID >= 0) return current_voxel_ID;
    return current_voxel_ID + num_sections;
}

[[nodiscard]] static std::vector<svr::SphericalVoxel>
svr::sphericalCoordinateVoxelTraversal(const Ray &ray, const SphericalVoxelGrid &grid, double t_begin,
                                       double t_end, double tol) noexcept {
    std::vector<svr::SphericalVoxel> voxels;
    voxels.reserve(grid.numRadialVoxels() * 3); // magic number

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
    const double discriminant = current_r * current_r - (ray_sphere_vector_dot - v * v);
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
    const double t_grid_end = std::max(t1, t2);

    /* TRAVERSAL PHASE */
    double t = t_begin;
    t_end = std::min(t_grid_end, t_end);
    bool previous_transition_flag = false;
    assert(false);

    return voxels;
}