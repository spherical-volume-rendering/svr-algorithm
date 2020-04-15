#include "../spherical_volume_rendering_util.h"

// Profiling over 4 cases:
// (1) Traversal from bottom left to upper right.
// (2) Traversal parallel to y-axis.
// (3) Traversal parallel to z-axis.
// (4) Traversal parallel to x-axis.
int main() {
    // Traversal from bottom left to upper right.
    const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
    const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
    const BoundVec3 sphere_center(0.0, 0.0, 0.0);
    const double sphere_max_radius = 10000.0;
    const std::size_t num_radial_sections = 10000;
    const std::size_t num_angular_sections = 10000;
    const std::size_t num_azimuthal_sections = 10000;
    const SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                  num_angular_sections,
                                  num_azimuthal_sections, sphere_center, sphere_max_radius);
    const BoundVec3 ray_origin(-13000.0, -13000.0, -13000.0);
    const FreeVec3 ray_direction(1.0, 1.0, 1.0);
    const Ray ray(ray_origin, ray_direction);
    const double t_begin = 0.0;
    const double t_end = 100000.0;
    const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);

    // Traversal parallel to the y-axis.
    const BoundVec3 min_bound2(-20000.0, -20000.0, -20000.0);
    const BoundVec3 max_bound2(20000.0, 20000.0, 20000.0);
    const BoundVec3 sphere_center2(0.0, 0.0, 0.0);
    const double sphere_max_radius2 = 10000.0;
    const std::size_t num_radial_sections2 = 10000;
    const std::size_t num_angular_sections2 = 10000;
    const std::size_t num_azimuthal_sections2 = 10000;
    const SphericalVoxelGrid grid2(min_bound2, max_bound2, num_radial_sections2,
                                  num_angular_sections2,
                                  num_azimuthal_sections2, sphere_center2, sphere_max_radius2);
    const BoundVec3 ray_origin2(10.0, -13000.0, 0.0);
    const FreeVec3 ray_direction2(0.0, 1.0, 0.0);
    const Ray ray2(ray_origin2, ray_direction2);
    const double t_begin2 = 0.0;
    const double t_end2 = 100000.0;
    const auto actual_voxels2 = sphericalCoordinateVoxelTraversal(ray2, grid2, t_begin2, t_end2);

    // Traversal parallel to the z-axis.
    const BoundVec3 min_bound3(-20000.0, -20000.0, -20000.0);
    const BoundVec3 max_bound3(20000.0, 20000.0, 20000.0);
    const BoundVec3 sphere_center3(0.0, 0.0, 0.0);
    const double sphere_max_radius3 = 10000.0;
    const std::size_t num_radial_sections3 = 10000;
    const std::size_t num_angular_sections3 = 10000;
    const std::size_t num_azimuthal_sections3 = 10000;
    const SphericalVoxelGrid grid3(min_bound3, max_bound3, num_radial_sections3,
                                   num_angular_sections3,
                                   num_azimuthal_sections3, sphere_center3, sphere_max_radius3);
    const BoundVec3 ray_origin3(10.0, 0.0, -13000.0);
    const FreeVec3 ray_direction3(0.0, 0.0, 1.0);
    const Ray ray3(ray_origin3, ray_direction3);
    const double t_begin3 = 0.0;
    const double t_end3 = 100000.0;
    const auto actual_voxels3 = sphericalCoordinateVoxelTraversal(ray3, grid3, t_begin3, t_end3);

    // Traversal parallel to the x-axis.
    const BoundVec3 min_bound4(-20000.0, -20000.0, -20000.0);
    const BoundVec3 max_bound4(20000.0, 20000.0, 20000.0);
    const BoundVec3 sphere_center4(0.0, 0.0, 0.0);
    const double sphere_max_radius4 = 10000.0;
    const std::size_t num_radial_sections4 = 10000;
    const std::size_t num_angular_sections4 = 10000;
    const std::size_t num_azimuthal_sections4 = 10000;
    const SphericalVoxelGrid grid4(min_bound4, max_bound4, num_radial_sections4,
                                  num_angular_sections4,
                                  num_azimuthal_sections4, sphere_center4, sphere_max_radius4);
    const BoundVec3 ray_origin4(-13000.0, 10.0, 10.0);
    const FreeVec3 ray_direction4(1.0, 0.0, 0.0);
    const Ray ray4(ray_origin4, ray_direction4);
    const double t_begin4 = 0.0;
    const double t_end4 = 100000.0;
    const auto actual_voxels4 = sphericalCoordinateVoxelTraversal(ray4, grid4, t_begin4, t_end4);
}

