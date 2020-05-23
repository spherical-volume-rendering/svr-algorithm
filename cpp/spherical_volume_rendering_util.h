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

  inline bool operator==(SphericalVoxel v) const noexcept {
    return this->radial == v.radial && this->polar == v.polar &&
           this->azimuthal == v.azimuthal;
  }
};

// A spherical coordinate voxel traversal algorithm. The algorithm traces the
// ray with unit direction over the spherical voxel grid provided. Returns a
// vector of the spherical coordinate voxels traversed. max_t is the unitized
// time for which the ray may travel. It is used in the following linear
// function: ray_travel_duration = time_to_entrance + sphere.diameter() * max_t
// If the ray origin is within the sphere, then time_to_entrance is 0. Its
// expected values are within bounds [0.0, 1.0]. For example, if max_t <= 0.0,
// then no voxels will be traversed. If max_t >= 1.0, then the entire sphere
// will be traversed.
// todo(cgyurgyik): Add description of sampler, data.
std::vector<SphericalVoxel> walkSphericalVolume(
    const Ray &ray, const svr::SphericalVoxelGrid &grid, double max_t,
    const std::function<void(const svr::SphericalVoxelGrid &grid,
                             const Ray &ray, double enter_t, double exit_t,
                             const svr::SphericalVoxel &sv, void *data)>
        &sampler = nullptr,
    void *data = nullptr) noexcept;

// Simplified parameters to Cythonize the function; implementation remains the
// same as above.
std::vector<SphericalVoxel> walkSphericalVolume(
    double *ray_origin, double *ray_direction, double *min_bound,
    double *max_bound, std::size_t num_radial_voxels,
    std::size_t num_polar_voxels, std::size_t num_azimuthal_voxels,
    double *sphere_center, double max_t) noexcept;

}  // namespace svr

#endif  // SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H