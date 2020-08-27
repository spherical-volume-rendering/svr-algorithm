#ifndef SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H
#define SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H

#include <vector>

#include "vec3.h"

namespace svr {

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

namespace {

constexpr double TAU = 2 * M_PI;

// Initializes the delta radii squared. These values are used for radial hit
// calculations in the main spherical volume algorithm. This calculates
// delta_radius^2 for num_radial_sections + 1 iterations. The delta radius value
// begins at max_radius, and subtracts delta_radius with each index.
// For example,
//
// Given: num_radial_voxels = 3, max_radius = 6, delta_radius = 2
// Returns: { 6*6, 4*4, 2*2, 0*0 }
std::vector<double> initializeDeltaRadiiSquared(
    const std::size_t num_radial_voxels, const double max_radius,
    const double delta_radius) noexcept {
  std::vector<double> delta_radii_squared(num_radial_voxels + 1);

  double current_delta_radius = max_radius;
  std::generate(delta_radii_squared.begin(), delta_radii_squared.end(),
                [&]() -> double {
                  const double old_delta_radius = current_delta_radius;
                  current_delta_radius -= delta_radius;
                  return old_delta_radius * old_delta_radius;
                });
  return delta_radii_squared;
}

// Returns a vector of TrigonometricValues for the given number of voxels.
// This begins with min_bound, and increments by a value of delta
// for num_voxels + 1 iterations. For example,
//
// Given: num_voxels = 2, min_bound = 0.0, delta = pi/2
// Returns: { {.cosine=1.0, .sine=0.0},
//            {.cosine=0.0, .sine=1.0},
//            {.cosine=1.0, .sine=0.0} }
std::vector<TrigonometricValues> initializeTrigonometricValues(
    const std::size_t num_voxels, const double min_bound,
    const double delta) noexcept {
  std::vector<TrigonometricValues> trig_values(num_voxels + 1);

  double radians = min_bound;
  std::generate(trig_values.begin(), trig_values.end(),
                [&]() -> TrigonometricValues {
                  const double cos = std::cos(radians);
                  const double sin = std::sin(radians);
                  radians += delta;
                  return {.cosine = cos, .sine = sin};
                });
  return trig_values;
}

// Returns a vector of maximum radius line segments for the given trigonometric
// values. The following predicate should be true:
// num_voxels + 1 == trig_values.size() + 1.
//
// The LineSegment points P1 and P2 are calculated with the following equation:
// .P1 = max_radius * trig_value.cosine + center.x().
// .P2 = max_radius * trig_value.sine + center.y().
std::vector<LineSegment> initializeMaxRadiusLineSegments(
    const std::size_t num_voxels, const BoundVec3 &center,
    const double max_radius,
    const std::vector<TrigonometricValues> &trig_values) noexcept {
  std::vector<LineSegment> line_segments(num_voxels + 1);
  std::transform(trig_values.cbegin(), trig_values.cend(),
                 line_segments.begin(),
                 [&](const TrigonometricValues &trig_value) -> LineSegment {
                   return {.P1 = max_radius * trig_value.cosine + center.x(),
                           .P2 = max_radius * trig_value.sine + center.y()};
                 });
  return line_segments;
}

// Initializes the vectors determined by the following calculation:
// sphere center - {X, Y, Z}, WHERE X, Y = P1, P2 for polar voxels.
std::vector<BoundVec3> initializeCenterToPolarPMaxVectors(
    const std::vector<LineSegment> &line_segments,
    const BoundVec3 &center) noexcept {
  std::vector<BoundVec3> center_to_pmax_vectors;
  center_to_pmax_vectors.reserve(line_segments.size());

  for (const auto &points : line_segments) {
    center_to_pmax_vectors.emplace_back(center -
                                        FreeVec3(points.P1, points.P2, 0.0));
  }
  return center_to_pmax_vectors;
}

// Similar to above, but uses:
// sphere center - {X, Y, Z}, WHERE X, Z = P1, P2 for azimuthal voxels.
std::vector<BoundVec3> initializeCenterToAzimuthalPMaxVectors(
    const std::vector<LineSegment> &line_segments,
    const BoundVec3 &center) noexcept {
  std::vector<BoundVec3> center_to_pmax_vectors;
  center_to_pmax_vectors.reserve(line_segments.size());

  for (const auto &points : line_segments) {
    center_to_pmax_vectors.emplace_back(center -
                                        FreeVec3(points.P1, 0.0, points.P2));
  }
  return center_to_pmax_vectors;
}

}  // namespace

// Represents a spherical voxel grid used for ray casting. The bounds of the
// grid are determined by min_bound and max_bound. The deltas are then
// determined by (max_bound.X - min_bound.X) / num_X_sections. To minimize
// calculation duplication, many calculations are completed once here and used
// each time a ray traverses the spherical voxel grid.
//
// Note that the grid system currently does not align with one would expect
// from spherical coordinates. We represent both polar and azimuthal within
// bounds [0, 2pi].
// TODO(cgyurgyik): Look into updating polar grid from [0, 2pi] -> [0, pi].
struct SphericalVoxelGrid {
 public:
  SphericalVoxelGrid(const SphereBound &min_bound, const SphereBound &max_bound,
                     std::size_t num_radial_sections,
                     std::size_t num_polar_sections,
                     std::size_t num_azimuthal_sections,
                     const BoundVec3 &sphere_center) noexcept
      : num_radial_sections_(num_radial_sections),
        num_polar_sections_(num_polar_sections),
        num_azimuthal_sections_(num_azimuthal_sections),
        sphere_center_(sphere_center),
        sphere_max_bound_polar_(max_bound.polar),
        sphere_min_bound_polar_(min_bound.polar),
        sphere_max_bound_azimuthal_(max_bound.azimuthal),
        sphere_min_bound_azimuthal_(min_bound.azimuthal),
        // TODO(cgyurgyik): Verify we want the sphere_max_radius to simply be
        // max_bound.radial.
        sphere_max_radius_(max_bound.radial),
        sphere_max_diameter_(sphere_max_radius_ * 2.0),
        delta_radius_((max_bound.radial - min_bound.radial) /
                      num_radial_sections),
        delta_theta_((max_bound.polar - min_bound.polar) / num_polar_sections),
        delta_phi_((max_bound.azimuthal - min_bound.azimuthal) /
                   num_azimuthal_sections),
        // TODO(cgyurgyik): Verify this is actually what we want for
        // 'max_radius'. The other option is simply using max_bound.radial
        delta_radii_sq_(initializeDeltaRadiiSquared(
            num_radial_sections,
            /*max_radius=*/max_bound.radial - min_bound.radial, delta_radius_)),
        polar_trig_values_(initializeTrigonometricValues(
            num_polar_sections, min_bound.polar, delta_theta_)),
        azimuthal_trig_values_(initializeTrigonometricValues(
            num_azimuthal_sections, min_bound.azimuthal, delta_phi_)),
        P_max_polar_(initializeMaxRadiusLineSegments(
            num_polar_sections, sphere_center, sphere_max_radius_,
            polar_trig_values_)),
        P_max_azimuthal_(initializeMaxRadiusLineSegments(
            num_azimuthal_sections, sphere_center, sphere_max_radius_,
            azimuthal_trig_values_)),
        center_to_polar_bound_vectors_(
            initializeCenterToPolarPMaxVectors(P_max_polar_, sphere_center)),
        center_to_azimuthal_bound_vectors_(
            initializeCenterToAzimuthalPMaxVectors(P_max_azimuthal_,
                                                   sphere_center)) {}

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

  // The delta radii squared calculated for use in radial hit calculations.
  const std::vector<double> delta_radii_sq_;

  // The trigonometric values calculated for the azimuthal and polar voxels.
  const std::vector<TrigonometricValues> azimuthal_trig_values_,
      polar_trig_values_;

  // The maximum radius line segments for polar and azimuthal voxels.
  const std::vector<LineSegment> P_max_polar_, P_max_azimuthal_;

  // The vectors represented by the vector sphere center - P_max[i].
  const std::vector<BoundVec3> center_to_polar_bound_vectors_,
      center_to_azimuthal_bound_vectors_;
};

}  // namespace svr

#endif  // SPHERICAL_VOLUME_RENDERING_SPHERICALVOXELGRID_H
