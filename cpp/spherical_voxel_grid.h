#ifndef SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H
#define SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H

#include <vector>

#include "vec3.h"

namespace svr {
constexpr double TAU = 2 * M_PI;

// Represents the boundary for the sphere. This is used to determine the minimum
// and maximum boundaries for a sectored traversal.
struct SphereBound {
  double radial;
  double polar;
  double azimuthal;
};

// Represents a line segment that is used for the points of intersections
// between the lines corresponding to voxel boundaries and a given radial voxel.
struct LineSegment {
  double P1;
  double P2;
};

// The trigonometric values for a given radian.
struct TrigonometricValues {
  double cosine;
  double sine;
};

// Represents a spherical voxel grid used for ray casting. The bounds of the
// grid are determined by min_bound and max_bound. The deltas are then
// determined by (max_bound.X - min_bound.X) / num_X_sections. To minimize
// calculation duplication, many calculations are completed once here and used
// each time a ray traverses the spherical voxel grid.
struct SphericalVoxelGrid {
 public:
  SphericalVoxelGrid(const SphereBound &min_bound, const SphereBound &max_bound,
                     std::size_t num_radial_sections,
                     std::size_t num_polar_sections,
                     std::size_t num_azimuthal_sections,
                     const BoundVec3 &sphere_center)
      : num_radial_sections_(num_radial_sections),
        num_polar_sections_(num_polar_sections),
        num_azimuthal_sections_(num_azimuthal_sections),
        sphere_center_(sphere_center),
        sphere_max_bound_polar_(max_bound.polar),
        sphere_min_bound_polar_(min_bound.polar),
        sphere_max_bound_azimuthal_(max_bound.azimuthal),
        sphere_min_bound_azimuthal_(min_bound.azimuthal),
        sphere_max_radius_(max_bound.radial),
        sphere_max_diameter_(sphere_max_radius_ * 2.0),
        delta_radius_((max_bound.radial - min_bound.radial) /
                      num_radial_sections),
        delta_theta_((max_bound.polar - min_bound.polar) / num_polar_sections),
        delta_phi_((max_bound.azimuthal - min_bound.azimuthal) /
                   num_azimuthal_sections) {
    double current_delta_radius = max_bound.radial - min_bound.radial;
    delta_radii_sq_.resize(num_radial_sections + 1);
    std::generate(delta_radii_sq_.begin(), delta_radii_sq_.end(),
                  [&]() -> double {
                    const double old_delta_radius = current_delta_radius;
                    current_delta_radius -= delta_radius_;
                    return old_delta_radius * old_delta_radius;
                  });
    P_max_polar_.resize(num_polar_sections + 1);
    P_max_azimuthal_.resize(num_azimuthal_sections + 1);
    center_to_polar_bound_vectors_.reserve(num_polar_sections + 1);
    center_to_azimuthal_bound_vectors_.reserve(num_azimuthal_sections + 1);
    if (num_polar_sections == num_azimuthal_sections) {
      double radians = 0.0;
      polar_trig_values_.resize(num_polar_sections + 1);
      std::generate(polar_trig_values_.begin(), polar_trig_values_.end(),
                    [&]() -> TrigonometricValues {
                      const double cos = std::cos(radians);
                      const double sin = std::sin(radians);
                      radians += delta_theta_;
                      return {.cosine = cos, .sine = sin};
                    });
      std::transform(polar_trig_values_.cbegin(), polar_trig_values_.cend(),
                     P_max_polar_.begin(), P_max_azimuthal_.begin(),
                     [&](const TrigonometricValues &tv,
                         LineSegment &ang_LS) -> LineSegment {
                       const double px_max_value =
                           sphere_max_radius_ * tv.cosine + sphere_center.x();
                       const double max_radius_times_s =
                           sphere_max_radius_ * tv.sine;
                       ang_LS = {.P1 = px_max_value,
                                 .P2 = max_radius_times_s + sphere_center.y()};
                       return {.P1 = px_max_value,
                               .P2 = max_radius_times_s + sphere_center.z()};
                     });
      for (const auto &points : P_max_polar_) {
        const BoundVec3 v =
            sphere_center - FreeVec3(points.P1, points.P2, points.P2);
        center_to_polar_bound_vectors_.push_back(v);
        center_to_azimuthal_bound_vectors_.push_back(v);
      }
      return;
    }
    double radians = 0.0;
    polar_trig_values_.resize(num_polar_sections + 1);
    std::generate(polar_trig_values_.begin(), polar_trig_values_.end(),
                  [&]() -> TrigonometricValues {
                    const double cos = std::cos(radians);
                    const double sin = std::sin(radians);
                    radians += delta_theta_;
                    return {.cosine = cos, .sine = sin};
                  });
    radians = 0.0;
    azimuthal_trig_values_.resize(num_azimuthal_sections + 1);
    std::generate(azimuthal_trig_values_.begin(), azimuthal_trig_values_.end(),
                  [&]() -> TrigonometricValues {
                    const double cos = std::cos(radians);
                    const double sin = std::sin(radians);
                    radians += delta_phi_;
                    return {.cosine = cos, .sine = sin};
                  });
    std::transform(
        polar_trig_values_.cbegin(), polar_trig_values_.cend(),
        P_max_polar_.begin(),
        [&](const TrigonometricValues &ang_tv) -> LineSegment {
          return {.P1 = sphere_max_radius_ * ang_tv.cosine + sphere_center.x(),
                  .P2 = sphere_max_radius_ * ang_tv.sine + sphere_center.y()};
        });
    std::transform(
        azimuthal_trig_values_.cbegin(), azimuthal_trig_values_.cend(),
        P_max_azimuthal_.begin(),
        [&](const TrigonometricValues &azi_tv) -> LineSegment {
          return {.P1 = sphere_max_radius_ * azi_tv.cosine + sphere_center.x(),
                  .P2 = sphere_max_radius_ * azi_tv.sine + sphere_center.z()};
        });
    for (const auto &points : P_max_polar_) {
      center_to_polar_bound_vectors_.emplace_back(
          sphere_center - FreeVec3(points.P1, points.P2, 0.0));
    }
    for (const auto &points : P_max_azimuthal_) {
      center_to_azimuthal_bound_vectors_.emplace_back(
          sphere_center - FreeVec3(points.P1, 0.0, points.P2));
    }
  }

  inline std::size_t numRadialSections() const noexcept {
    return this->num_radial_sections_;
  }

  inline std::size_t numPolarSections() const noexcept {
    return this->num_polar_sections_;
  }

  inline std::size_t numAzimuthalSections() const noexcept {
    return this->num_azimuthal_sections_;
  }

  inline double sphereMaxBoundPolar() const noexcept {
    return this->sphere_max_bound_polar_;
  }

  inline double sphereMinBoundPolar() const noexcept {
    return this->sphere_min_bound_polar_;
  }

  inline double sphereMaxBoundAzi() const noexcept {
    return this->sphere_max_bound_azimuthal_;
  }

  inline double sphereMinBoundAzi() const noexcept {
    return this->sphere_min_bound_azimuthal_;
  }

  inline double sphereMaxRadius() const noexcept {
    return this->sphere_max_radius_;
  }

  inline double sphereMaxDiameter() const noexcept {
    return this->sphere_max_diameter_;
  }

  inline const BoundVec3 &sphereCenter() const noexcept {
    return this->sphere_center_;
  }

  inline double deltaRadius() const noexcept { return delta_radius_; }

  inline double deltaPhi() const noexcept { return delta_phi_; }

  inline double deltaTheta() const noexcept { return delta_theta_; }

  inline double deltaRadiiSquared(std::size_t i) const noexcept {
    return this->delta_radii_sq_[i];
  }

  inline const LineSegment &pMaxPolar(std::size_t i) const noexcept {
    return this->P_max_polar_[i];
  }

  inline const std::vector<LineSegment> &pMaxPolar() const noexcept {
    return this->P_max_polar_;
  }

  inline const BoundVec3 &centerToPolarBound(std::size_t i) const noexcept {
    return this->center_to_polar_bound_vectors_[i];
  }

  inline const LineSegment &pMaxAzimuthal(std::size_t i) const noexcept {
    return this->P_max_azimuthal_[i];
  }

  inline const std::vector<LineSegment> &pMaxAzimuthal() const noexcept {
    return this->P_max_azimuthal_;
  }

  inline const BoundVec3 &centerToAzimuthalBound(std::size_t i) const noexcept {
    return this->center_to_azimuthal_bound_vectors_[i];
  }

  inline const std::vector<TrigonometricValues> &polarTrigValues()
      const noexcept {
    return polar_trig_values_;
  }

  inline const std::vector<TrigonometricValues> &azimuthalTrigValues()
      const noexcept {
    return azimuthal_trig_values_;
  }

 private:
  // The number of radial, polar, and azimuthal voxels.
  const std::size_t num_radial_sections_, num_polar_sections_,
      num_azimuthal_sections_;

  // The center of the sphere.
  const BoundVec3 sphere_center_;

  // The maximum polar bound of the sphere.
  const double sphere_max_bound_polar_;

  // The minimum polar bound of the sphere.
  const double sphere_min_bound_polar_;

  // The maximum azimuthal bound of the sphere.
  const double sphere_max_bound_azimuthal_;

  // The minimum azimuthal bound of the sphere.
  const double sphere_min_bound_azimuthal_;

  // The maximum radius of the sphere.
  const double sphere_max_radius_;

  // The maximum diamater of the sphere.
  const double sphere_max_diameter_;

  // The maximum sphere radius divided by the number of radial sections.
  const double delta_radius_;

  // 2 * PI divided by X, where X is the number of polar and number of azimuthal
  // sections respectively.
  const double delta_theta_, delta_phi_;

  // The delta radii squared ranging from 0...num_radial_voxels.
  std::vector<double> delta_radii_sq_;

  // The maximum radius line segments for polar voxels.
  std::vector<LineSegment> P_max_polar_;

  // The maximum radius line segments for azimuthal voxels.
  std::vector<LineSegment> P_max_azimuthal_;

  // The trigonometric values for each delta theta.
  std::vector<TrigonometricValues> polar_trig_values_;

  // The trigonometric values for each delta phi. In the case where delta theta
  // is equal to delta phi, this is left uninitialized and polar_trig_values_ is
  // used.
  std::vector<TrigonometricValues> azimuthal_trig_values_;

  // The vectors represented by sphere_center_ - P_max_polar_[i].
  std::vector<BoundVec3> center_to_polar_bound_vectors_;

  // The vectors represented by sphere_center_ - P_max_azimuthal_[i].
  std::vector<BoundVec3> center_to_azimuthal_bound_vectors_;
};

}  // namespace svr

#endif  // SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H
