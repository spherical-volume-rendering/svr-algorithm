#ifndef SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H
#define SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H

#include "ray.h"
#include "vec3.h"
#include "spherical_voxel_grid.h"
#include <vector>

namespace svr {

    // Represents a spherical voxel coordinate.
    struct SphericalVoxel {
        int radial_voxel;
        int angular_voxel;
        int azimuthal_voxel;
    };

    // Represents a bound for the sphere using voxel identifications.
    struct VoxelBound {
        int radial_voxel;
        int angular_voxel;
        int azimuthal_voxel;
    };

    // A spherical coordinate voxel traversal algorithm. The algorithm traces the ray over the spherical voxel grid
    // provided. t_begin is the time the ray begins, and t_end is the time at which the ray ends. Returns a
    // vector of the spherical coordinate voxels traversed. Does not assume that the ray direction is normalized.
    // Recall that if, for example, a radial hit occurs, the azimuthal and angular voxels will remain the
    // same as before. This applies for each traversal type. To determine the boundaries of intersection,
    // VoxelBound min_bound & max_bound are used. Examples of boundary use are provided below.
    //
    // Example 1: If one wants to traverse the entire sphere without limitation:
    //            VoxelBound min_bound = {.radial_voxel = 1, .angular_voxel = 0, .azimuthal_voxel = 0};
    //            VoxelBound max_bound = {.radial_voxel = num_radial_sections,
    //                                    .angular_voxel = num_angular_sections - 1,
    //                                    .azimuthal_voxel = num_azimuthal_sections - 1};
    //
    // Example 2: If one wants to traverse a sector of the sphere, say the first quadrant:
    //            If num_radial_sections = num_angular_sections = num_azimuthal_sections = 4, then
    //            VoxelBound min_bound = {.radial_voxel = 1, .angular_voxel = 0, .azimuthal_voxel = 0};
    //            VoxelBound max_bound = {.radial_voxel = 3, .angular_voxel = 0, .azimuthal_voxel = 0};
    //
    // Note: We assume for each voxel v, min_bound.v < max_bound.v. Therefore, if one wanted to traverse the section
    //       of a circle representing (pi/4, -pi/4) -> (pi/4, pi/4), this would need to be broken into two separate
    //       function calls. Similarly, for a voxel with 8 angular sections, if one wants to traverse from
    //       (angular section 7) -> (angular section 0), this would need to be broken into two separate function calls.
    std::vector<SphericalVoxel> walkSphericalVolume(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                                    const svr::VoxelBound &min_bound,
                                                    const svr::VoxelBound &max_bound,
                                                    double t_begin, double t_end) noexcept;

    // Simplified parameters to Cythonize the function; implementation remains the same as above.
    std::vector<SphericalVoxel> walkSphericalVolume(double *ray_origin, double *ray_direction,
                                                    int *voxel_min_bound, int *voxel_max_bound,
                                                    std::size_t num_radial_voxels, std::size_t num_angular_voxels,
                                                    std::size_t num_azimuthal_voxels, double *sphere_center,
                                                    double sphere_max_radius, double t_begin, double t_end) noexcept;

} // namespace svr

#endif //SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H
