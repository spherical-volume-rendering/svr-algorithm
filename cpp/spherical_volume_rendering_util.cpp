#include "spherical_volume_rendering_util.h"
#include <vector>

// Takes the absolute value of 'value' and determines if it is less than the tolerance.
// In this case, we assume tolerance is an arbitrary small value.
// This is used when floating point mathematics carries some rounding errors.
inline bool isNearZero(double value, double tolerance) {
    return std::abs(value) < tolerance;
}

[[nodiscard]] static std::vector<svr::SphericalVoxel>
svr::sphericalCoordinateVoxelTraversal(const Ray &ray, const SphericalVoxelGrid &grid, double t_begin,
                                       double t_end, double tol) noexcept {
    std::vector<svr::SphericalVoxel> voxels;
    voxels.reserve(grid.numRadialVoxels() * 3);

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

    // TODO(cgyurgyik): Finish voxel ID theta calculation.
    assert(false); // UNFINISHED.

    return voxels;
}