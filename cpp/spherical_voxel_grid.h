#ifndef SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H
#define SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H

#include "vec3.h"

// Represents a 3-dimensional spherical voxel grid. The minimum and maximum bounds [min_bound_, max_bound_]
// contain the entirety of the sphere that is to be traversed.
// Requires:
//   max_bound > min_bound
//   num_radial_voxels > 0
//   num_angular_voxels > 0
//   num_azimuthal_voxels > 0
//   sphere_max_radius > 0.0
struct SphericalVoxelGrid {
public:
    constexpr SphericalVoxelGrid(const BoundVec3& min_bound, const BoundVec3& max_bound,
            size_t num_radial_voxels, size_t num_angular_voxels, size_t num_azimuthal_voxels,
            const BoundVec3& sphere_center, float sphere_max_radius) :
            min_bound_{min_bound},
            max_bound_{max_bound},
            num_radial_voxels_{num_radial_voxels},
            num_angular_voxels_{num_angular_voxels},
            num_azimuthal_voxels_{num_angular_voxels},
            sphere_center_{sphere_center},
            sphere_max_radius_{sphere_max_radius},
            inv_num_radial_voxels_{1.0 / num_radial_voxels},
            inv_num_angular_voxels_{1.0 / num_angular_voxels},
            inv_num_azimuthal_voxels_{1.0 / num_azimuthal_voxels},
            delta_radius_{sphere_max_radius * inv_num_radial_voxels_},
            delta_theta_{2 * M_PI * inv_num_angular_voxels_},
            delta_phi_{2 * M_PI * inv_num_azimuthal_voxels_}
            {}

    [[nodiscard]] inline constexpr size_t numRadialVoxels() const { return num_radial_voxels_; }

    [[nodiscard]] inline constexpr size_t numAngularVoxels() const { return num_angular_voxels_; }

    [[nodiscard]] inline constexpr size_t numAzimuthalVoxels() const { return num_azimuthal_voxels_; }

    [[nodiscard]] inline constexpr double invNumRadialVoxels() const { return inv_num_radial_voxels_; }

    [[nodiscard]] inline constexpr double invNumAngularVoxels() const { return inv_num_angular_voxels_; }

    [[nodiscard]] inline constexpr double invNumAzimuthalVoxels() const { return inv_num_azimuthal_voxels_; }

    [[nodiscard]] inline constexpr BoundVec3 minBound() const { return min_bound_; }

    [[nodiscard]] inline constexpr BoundVec3 maxBound() const { return max_bound_; }

    [[nodiscard]] inline constexpr double sphereMaxRadius() const { return sphere_max_radius_; }

    [[nodiscard]] inline constexpr BoundVec3 sphereCenter() const { return sphere_center_; }

    [[nodiscard]] inline constexpr double deltaRadius() const { return delta_radius_; }

    [[nodiscard]] inline constexpr double deltaTheta() const { return delta_theta_; }

    [[nodiscard]] inline constexpr double deltaPhi() const { return delta_phi_; }

private:
    // The minimum bound vector of the voxel grid.
    const BoundVec3 min_bound_;

    // The maximum bound vector of the voxel grid.
    const BoundVec3 max_bound_;

    // The number of radial, angular, and azimuthal voxels.
    const size_t num_radial_voxels_, num_angular_voxels_, num_azimuthal_voxels_;

    // Similar to num_x_voxels, but 1 / x, where x is the number of voxels.
    const double inv_num_radial_voxels_, inv_num_angular_voxels_, inv_num_azimuthal_voxels_;

    // The center of the sphere.
    const BoundVec3 sphere_center_;

    // The maximum radius of the sphere.
    const double sphere_max_radius_;

    // The maximum sphere radius divided by the number of radial sections.
    const double delta_radius_;

    // Sphere's area divided by the number of angular sections.
    const double delta_theta_;

    // Sphere's area divided by the number of azimuthal sections;
    const double delta_phi_;
};

#endif //SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H
