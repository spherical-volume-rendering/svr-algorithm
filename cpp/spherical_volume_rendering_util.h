#ifndef SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H
#define SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H

#include <vector>

#include "ray.h"
#include "spherical_voxel_grid.h"
#include "vec3.h"

namespace svr {

// Represents a spherical voxel coordinate.
struct SphericalVoxel {
  int radial;
  int polar;
  int azimuthal;
};

// A spherical coordinate voxel traversal algorithm. The algorithm traces the
// ray over the spherical voxel grid provided. max_t is the maximum unitized time at
// until which the ray traverses. Returns a vector of the spherical coordinate
// voxels traversed.
// todo: Explain max_t.
std::vector<SphericalVoxel> walkSphericalVolume(
    const Ray &ray, const svr::SphericalVoxelGrid &grid, double max_t) noexcept;

// Simplified parameters to Cythonize the function; implementation remains the
// same as above.
std::vector<SphericalVoxel> walkSphericalVolume(
    double *ray_origin, double *ray_direction, double *min_bound,
    double *max_bound, std::size_t num_radial_voxels,
    std::size_t num_polar_voxels, std::size_t num_azimuthal_voxels,
    double *sphere_center, double max_t) noexcept;

}  // namespace svr

#endif  // SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H