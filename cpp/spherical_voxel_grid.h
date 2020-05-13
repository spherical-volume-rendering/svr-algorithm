#ifndef SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H
#define SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H

#include "vec3.h"
#include <vector>

namespace svr {
    constexpr double TAU = 2 * M_PI;

    // Represents the boundary for the sphere. This is used to determine the minimum and maximum
    // boundaries for a sectored traversal.
    struct SphereBound {
        double radial;
        double polar;
        double azimuthal;
    };

    // Represents a line segment that is used for the points of intersections between
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


    // Represents a spherical voxel grid used for ray casting. The bounds of the grid are determined by min_bound
    // and max_bound. The deltas are then determined by (max_bound.X - min_bound.X) / num_X_sections. To minimize
    // calculation duplication, many calculations are completed once here and used each time a ray traverses the
    // spherical voxel grid.
    struct SphericalVoxelGrid {
    public:
        SphericalVoxelGrid(const SphereBound &min_bound, const SphereBound &max_bound,
                std::size_t num_radial_voxels, std::size_t num_polar_voxels, std::size_t num_azimuthal_voxels,
                const BoundVec3 &sphere_center) :
                num_radial_voxels_(num_radial_voxels),
                num_polar_voxels_(num_polar_voxels),
                num_azimuthal_voxels_(num_azimuthal_voxels),
                sphere_center_(sphere_center),
                sphere_max_radius_(max_bound.radial),
                delta_radius_((max_bound.radial - min_bound.radial) / num_radial_voxels),
                delta_theta_((max_bound.polar - min_bound.polar) / num_polar_voxels),
                delta_phi_((max_bound.azimuthal - min_bound.azimuthal) / num_azimuthal_voxels),
                inverse_delta_radius_(1.0 / delta_radius_),
                inverse_delta_theta_(1.0 / delta_theta_),
                inverse_delta_phi_(1.0 / delta_phi_) {

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

            P_max_polar_.resize(num_polar_voxels + 1);
            P_max_azimuthal_.resize(num_azimuthal_voxels + 1);
            center_to_polar_bound_vectors_.reserve(num_polar_voxels + 1);
            center_to_azimuthal_bound_vectors_.reserve(num_azimuthal_voxels + 1);
            if (num_polar_voxels == num_azimuthal_voxels) {
                double radians = 0.0;
                polar_trig_values_.resize(num_polar_voxels + 1);
                std::generate(polar_trig_values_.begin(), polar_trig_values_.end(),
                              [&]() -> TrigonometricValues {
                                  const double cos = std::cos(radians);
                                  const double sin = std::sin(radians);
                                  radians += delta_theta_;
                                  return {.cosine=cos, .sine=sin};
                              });
                std::transform(polar_trig_values_.cbegin(), polar_trig_values_.cend(),
                               P_max_polar_.begin(),
                               P_max_azimuthal_.begin(),
                               [&](const TrigonometricValues &tv, LineSegment &ang_LS) -> LineSegment {
                                   const double px_max_value = sphere_max_radius_ * tv.cosine + sphere_center.x();
                                   const double max_radius_times_s = sphere_max_radius_ * tv.sine;
                                   ang_LS = {.P1=px_max_value, .P2=max_radius_times_s + sphere_center.y()};
                                   return {.P1=px_max_value, .P2=max_radius_times_s + sphere_center.z()};
                               });
                for (const auto& points : P_max_polar_) {
                    const BoundVec3 v = sphere_center - FreeVec3(points.P1, points.P2, points.P2);
                    center_to_polar_bound_vectors_.push_back(v);
                    center_to_azimuthal_bound_vectors_.push_back(v);
                }
                return;
            }

            double radians = 0.0;
            polar_trig_values_.resize(num_polar_voxels + 1);
            std::generate(polar_trig_values_.begin(), polar_trig_values_.end(), [&]() -> TrigonometricValues {
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
            std::transform(polar_trig_values_.cbegin(), polar_trig_values_.cend(), P_max_polar_.begin(),
                           [&](const TrigonometricValues &ang_tv) -> LineSegment {
                               return {.P1=sphere_max_radius_ * ang_tv.cosine + sphere_center.x(),
                                       .P2=sphere_max_radius_ * ang_tv.sine + sphere_center.y()};
                           });
            std::transform(azimuthal_trig_values_.cbegin(), azimuthal_trig_values_.cend(), P_max_azimuthal_.begin(),
                           [&](const TrigonometricValues &azi_tv) -> LineSegment {
                               return {.P1=sphere_max_radius_ * azi_tv.cosine + sphere_center.x(),
                                       .P2=sphere_max_radius_ * azi_tv.sine + sphere_center.z()};
                           });
            for (const auto& points : P_max_polar_) {
                center_to_polar_bound_vectors_.emplace_back(sphere_center - FreeVec3(points.P1, points.P2, 0.0));
            }
            for (const auto& points : P_max_azimuthal_) {
                center_to_azimuthal_bound_vectors_.emplace_back(sphere_center - FreeVec3(points.P1, 0.0, points.P2));
            }
        }

        inline std::size_t numRadialVoxels() const noexcept { return this->num_radial_voxels_; }

        inline std::size_t numPolarVoxels() const noexcept { return this->num_polar_voxels_; }

        inline std::size_t numAzimuthalVoxels() const noexcept { return this->num_azimuthal_voxels_; }

        inline double sphereMaxRadius() const noexcept { return this->sphere_max_radius_; }

        inline const BoundVec3 &sphereCenter() const noexcept { return this->sphere_center_; }

        inline double deltaRadii(std::size_t i) const noexcept { return this->delta_radii_[i]; }

        inline double deltaRadiiSquared(std::size_t i) const noexcept { return this->delta_radii_sq_[i]; }

        inline const std::vector<double> &deltaRadii() const noexcept { return this->delta_radii_; }

        inline const std::vector<double> &deltaRadiiSquared() const noexcept { return this->delta_radii_sq_; }

        inline const LineSegment &pMaxPolar(std::size_t i) const noexcept { return this->P_max_polar_[i]; }

        inline const std::vector<LineSegment> &pMaxPolar() const noexcept { return this->P_max_polar_; }

        inline const BoundVec3 &
        centerToPolarBound(std::size_t i) const noexcept { return this->center_to_polar_bound_vectors_[i]; }

        inline const LineSegment &pMaxAzimuthal(std::size_t i) const noexcept { return this->P_max_azimuthal_[i]; }

        inline const std::vector<LineSegment> &pMaxAzimuthal() const noexcept { return this->P_max_azimuthal_; }

        inline const BoundVec3 &
        centerToAzimuthalBound(std::size_t i) const noexcept { return this->center_to_azimuthal_bound_vectors_[i]; }

        inline const std::vector<TrigonometricValues> &polarTrigValues() const noexcept { return polar_trig_values_; }

        inline const std::vector<TrigonometricValues> &
        azimuthalTrigValues() const noexcept { return azimuthal_trig_values_; }

    private:
        // The number of radial, polar, and azimuthal voxels.
        const std::size_t num_radial_voxels_, num_polar_voxels_, num_azimuthal_voxels_;

        // The center of the sphere.
        const BoundVec3 sphere_center_;

        // The maximum radius of the sphere.
        const double sphere_max_radius_;

        // The maximum sphere radius divided by the number of radial sections.
        const double delta_radius_;

        // 2 * PI divided by X, where X is the number of polar and number of azimuthal sections respectively.
        const double delta_theta_, delta_phi_;

        // The inverse of the deltas.
        const double inverse_delta_radius_, inverse_delta_theta_, inverse_delta_phi_;

        // The delta radii ranging from 0...num_radial_voxels.
        std::vector<double> delta_radii_;

        // The delta radii squared ranging from 0...num_radial_voxels.
        std::vector<double> delta_radii_sq_;

        // The maximum radius line segments for polar voxels.
        std::vector<LineSegment> P_max_polar_;

        // The maximum radius line segments for azimuthal voxels.
        std::vector<LineSegment> P_max_azimuthal_;

        // The trigonometric values for each delta theta.
        std::vector<TrigonometricValues> polar_trig_values_;

        // The trigonometric values for each delta phi. In the case where delta theta is equal to delta phi,
        // this is left uninitialized and polar_trig_values_ is used.
        std::vector<TrigonometricValues> azimuthal_trig_values_;

        // The vectors represented by sphere_center_ - P_max_polar_[i].
        std::vector<BoundVec3> center_to_polar_bound_vectors_;

        // The vectors represented by sphere_center_ - P_max_azimuthal_[i].
        std::vector<BoundVec3> center_to_azimuthal_bound_vectors_;
    };

} // namespace svr

#endif //SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H