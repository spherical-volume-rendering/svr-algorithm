#ifndef SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H
#define SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H

#include "ray.h"
#include "vec3.h"
#include "spherical_voxel_grid.h"
#include <vector>

// Represents a spherical voxel coordinate.
    struct SphericalVoxel {
        size_t radial_voxel;
        size_t angular_voxel;
        size_t azimuthal_voxel;
    };

// A spherical coordinate voxel traversal algorithm. The algorithm traces the 'ray' over the spherical voxel grid
// provided. 't_begin' is the time the ray begins, and 't_end' is the time the ray ends.
// 'tol' is the tolerance allowed for error that carries with any floating point operations.
// Requires:
//    'ray' is a valid Ray.
//    'grid' is a valid SphericalVoxelGrid.
//    t_end > t_begin >= 0.0
//    'tol' >= 0.0
// Returns:
//    A vector of the spherical coordinate voxels traversed. Recall that if a radial hit occurs,
//    The azimuthal and angular voxels will remain the same as before.
std::vector<SphericalVoxel> sphericalCoordinateVoxelTraversal(const Ray &ray, const SphericalVoxelGrid &grid,
                                                              double t_begin, double t_end, double tol) noexcept;

// Simplified version to Cythonize.
std::vector<SphericalVoxel> sphericalCoordinateVoxelTraversalCy(double* ray_origin, double* ray_direction,
                                                              double* min_bound, double* max_bound,
                                                              size_t num_radial_voxels, size_t num_angular_voxels,
                                                              size_t num_azimuthal_voxels, double* sphere_center,
                                                              double sphere_max_radius, double t_begin, double t_end,
                                                              double tol) noexcept;

#endif //SPHERICAL_VOLUME_RENDERING_SPHERICALVOLUMERENDERINGUTIL_H
