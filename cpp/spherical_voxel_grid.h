#ifndef SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H
#define SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H

#include "vec3.h"
#include <vector>

namespace svr {
    constexpr double TAU = 2 * M_PI;

    // Represents a bound for the sphere using voxel identifications.
    struct VoxelBound {
        std::size_t radial_voxel;
        std::size_t angular_voxel;
        std::size_t azimuthal_voxel;
    };

    // Represents a line segment. This is used to represent the points of intersections between
    // the lines corresponding to voxel boundaries and a given radial voxel.
    struct LineSegment {
        double P1;
        double P2;
    };

    // The trigonometric values for a given radian.
    struct TrigonometricValues {
        double cosine;
        double sine;
    };

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
        SphericalVoxelGrid(std::size_t num_radial_voxels, std::size_t num_angular_voxels,
                std::size_t num_azimuthal_voxels, const BoundVec3 &sphere_center, double sphere_max_radius) :
                num_radial_voxels_(num_radial_voxels),
                num_angular_voxels_(num_angular_voxels),
                num_azimuthal_voxels_(num_azimuthal_voxels),
                sphere_center_(sphere_center),
                sphere_max_radius_(sphere_max_radius),
                delta_radius_(sphere_max_radius / num_radial_voxels),
                delta_theta_(TAU / num_angular_voxels),
                delta_phi_(TAU / num_azimuthal_voxels) {

            delta_radii_.resize(num_radial_voxels + 1);
            double current_delta_radius = delta_radius_ * num_radial_voxels;
            std::generate(delta_radii_.begin(), delta_radii_.end(),
                          [&]() -> double {
                              const double old_delta_radius = current_delta_radius;
                              current_delta_radius -= delta_radius_;
                              return old_delta_radius;
                          });
            delta_radii_sq_.resize(num_radial_voxels + 1);
            std::transform(delta_radii_.cbegin(), delta_radii_.cend(), delta_radii_sq_.begin(),
                           [](double dR) -> double { return dR * dR; });

            P_max_angular_.resize(num_angular_voxels + 1);
            P_max_azimuthal_.resize(num_azimuthal_voxels + 1);
            if (num_angular_voxels == num_azimuthal_voxels) {
                double radians = 0.0;
                angular_trig_values_.resize(num_angular_voxels + 1);
                std::generate(angular_trig_values_.begin(), angular_trig_values_.end(),
                              [&]() -> TrigonometricValues {
                                  const double cos = std::cos(radians);
                                  const double sin = std::sin(radians);
                                  radians += delta_theta_;
                                  return {.cosine=cos, .sine=sin};
                              });
                std::transform(angular_trig_values_.cbegin(), angular_trig_values_.cend(),
                               P_max_angular_.begin(),
                               P_max_azimuthal_.begin(),
                               [&](const TrigonometricValues &tv, LineSegment &ang_LS) -> LineSegment {
                                   const double px_max_value = sphere_max_radius * tv.cosine + sphere_center.x();
                                   const double max_radius_times_s = sphere_max_radius * tv.sine;
                                   ang_LS = {.P1=px_max_value, .P2=max_radius_times_s + sphere_center.y()};
                                   return {.P1=px_max_value, .P2=max_radius_times_s + sphere_center.z()};
                               });
                return;
            }

            double radians = 0.0;
            angular_trig_values_.resize(num_angular_voxels + 1);
            std::generate(angular_trig_values_.begin(), angular_trig_values_.end(), [&]() -> TrigonometricValues {
                const double cos = std::cos(radians);
                const double sin = std::sin(radians);
                radians += delta_theta_;
                return {.cosine=cos, .sine=sin};
            });
            radians = 0.0;
            azimuthal_trig_values_.resize(num_azimuthal_voxels + 1);
            std::generate(azimuthal_trig_values_.begin(), azimuthal_trig_values_.end(), [&]() -> TrigonometricValues {
                const double cos = std::cos(radians);
                const double sin = std::sin(radians);
                radians += delta_phi_;
                return {.cosine=cos, .sine=sin};
            });
            std::transform(angular_trig_values_.cbegin(), angular_trig_values_.cend(), P_max_angular_.begin(),
                           [&](const TrigonometricValues &ang_tv) -> LineSegment {
                               return {.P1=sphere_max_radius * ang_tv.cosine + sphere_center.x(),
                                       .P2=sphere_max_radius * ang_tv.sine + sphere_center.y()};
                           });
            std::transform(azimuthal_trig_values_.cbegin(), azimuthal_trig_values_.cend(), P_max_azimuthal_.begin(),
                           [&](const TrigonometricValues &azi_tv) -> LineSegment {
                               return {.P1=sphere_max_radius * azi_tv.cosine + sphere_center.x(),
                                       .P2=sphere_max_radius * azi_tv.sine + sphere_center.z()};
                           });
        }

        inline std::size_t numRadialVoxels() const noexcept { return this->num_radial_voxels_; }

        inline std::size_t numAngularVoxels() const noexcept { return this->num_angular_voxels_; }

        inline std::size_t numAzimuthalVoxels() const noexcept { return this->num_azimuthal_voxels_; }

        inline double sphereMaxRadius() const noexcept { return this->sphere_max_radius_; }

        inline const BoundVec3 &sphereCenter() const noexcept { return this->sphere_center_; }

        inline double deltaRadius() const noexcept { return this->delta_radius_; }

        inline double deltaRadii(std::size_t i) const noexcept { return this->delta_radii_[i]; }

        inline double deltaRadiiSquared(std::size_t i) const noexcept { return this->delta_radii_sq_[i]; }

        inline const std::vector<double> &deltaRadii() const noexcept { return this->delta_radii_; }

        inline const std::vector<double> &deltaRadiiSquared() const noexcept { return this->delta_radii_sq_; }

        inline const LineSegment &pMaxAngular(std::size_t i) const noexcept { return this->P_max_angular_[i]; }

        inline const std::vector<LineSegment> &pMaxAngular() const noexcept { return this->P_max_angular_; }

        inline const LineSegment &pMaxAzimuthal(std::size_t i) const noexcept { return this->P_max_azimuthal_[i]; }

        inline const std::vector<LineSegment> &pMaxAzimuthal() const noexcept { return this->P_max_azimuthal_; }

        inline const std::vector<TrigonometricValues> &
        angularTrigValues() const noexcept { return angular_trig_values_; }

        inline const std::vector<TrigonometricValues> &
        azimuthalTrigValues() const noexcept { return azimuthal_trig_values_; }

    private:
        // The number of radial, angular, and azimuthal voxels.
        const std::size_t num_radial_voxels_, num_angular_voxels_, num_azimuthal_voxels_;

        // The center of the sphere.
        const BoundVec3 sphere_center_;

        // The maximum radius of the sphere.
        const double sphere_max_radius_;

        // The maximum sphere radius divided by the number of radial sections.
        const double delta_radius_;

        // 2 * PI divided by X, where X is the number of angular and number of azimuthal sections respectively.
        const double delta_theta_, delta_phi_;

        // The delta radii ranging from 0...num_radial_voxels.
        std::vector<double> delta_radii_;

        // The delta radii squared ranging from 0...num_radial_voxels.
        std::vector<double> delta_radii_sq_;

        // The maximum radius line segments for angular voxels.
        std::vector<LineSegment> P_max_angular_;

        // The maximum radius line segments for azimuthal voxels.
        std::vector<LineSegment> P_max_azimuthal_;

        // The trigonometric values for each delta theta.
        std::vector<TrigonometricValues> angular_trig_values_;

        // The trigonometric values for each delta phi. In the case where delta theta is equal to delta phi,
        // this is left uninitialized and angular_trig_values_ is used.
        std::vector<TrigonometricValues> azimuthal_trig_values_;
    };

} // namespace svr

#endif //SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H