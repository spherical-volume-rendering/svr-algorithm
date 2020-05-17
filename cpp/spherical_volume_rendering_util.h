#ifndef SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H
#define SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H

#include "ray.h"
#include "vec3.h"
#include "spherical_voxel_grid.h"
#include <vector>

namespace svr {

    // Represents a spherical voxel coordinate.
    struct SphericalVoxel {
        int radial;
        int polar;
        int azimuthal;
    };

    // A spherical coordinate voxel traversal algorithm. The algorithm traces the ray over the spherical voxel grid
    // provided. t_begin is the time the ray begins, and t_end is the time at which the ray ends. It does not
    // assumes the ray direction is normalized. Returns a vector of the spherical coordinate voxels traversed.
    std::vector<SphericalVoxel> walkSphericalVolume(const Ray &ray, const svr::SphericalVoxelGrid &grid,
                                                    double t_begin, double t_end) noexcept;

    // Simplified parameters to Cythonize the function; implementation remains the same as above.
    std::vector<SphericalVoxel> walkSphericalVolume(double *ray_origin, double *ray_direction,
                                                    double *min_bound, double *max_bound,
                                                    std::size_t num_radial_voxels,
                                                    std::size_t num_polar_voxels,
                                                    std::size_t num_azimuthal_voxels,
                                                    double *sphere_center, double t_begin, double t_end) noexcept;

} // namespace svr

#endif //SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H