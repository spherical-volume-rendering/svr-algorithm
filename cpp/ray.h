#ifndef SPHERICAL_VOLUME_RENDERING_RAY_H
#define SPHERICAL_VOLUME_RENDERING_RAY_H

#include "vec3.h"

// The indices for Vec3. For example, Vec3[0] returns the x-direction.
enum NonZeroDirectionIndex {
  X_DIRECTION = 0, Y_DIRECTION = 1, Z_DIRECTION = 2
};

// This struct provides parameters for division by a non-zero direction.
// This avoids the need to check if the direction is non-zero each time a function below is called.
struct NonZeroDirectionParameters {
  double inv_direction;
  double unit_direction;
  double origin;
  NonZeroDirectionIndex index;
};

// Encapsulates the functionality of a ray. This consists of two components, the origin of the ray, and the
// direction of the ray. To avoid checking for a non-zero direction upon each function call, these parameters are
// initialized upon construction.
struct Ray final {
    Ray(const BoundVec3 &origin, const FreeVec3 &direction)
            : origin_(origin), direction_(direction), unit_direction_(direction),
              inverse_direction_(FreeVec3(1.0 / direction.x(), 1.0 / direction.y(), 1.0 / direction.z())) {

      if (std::abs(direction.x()) > 0) {
        NZD_params_.origin = origin.x();
        NZD_params_.inv_direction = inverse_direction_.x();
        NZD_params_.unit_direction = unit_direction_.x();
        NZD_params_.index = X_DIRECTION;
        return;
      }
      if (std::abs(direction.y()) > 0) {
        NZD_params_.origin = origin.y();
        NZD_params_.inv_direction = inverse_direction_.y();
        NZD_params_.unit_direction = unit_direction_.y();
        NZD_params_.index = Y_DIRECTION;
        return;
      }
      NZD_params_.origin = origin.z();
      NZD_params_.inv_direction = inverse_direction_.z();
      NZD_params_.unit_direction = unit_direction_.z();
      NZD_params_.index = Z_DIRECTION;
    }

    // Represents the function p(t) = origin + t * direction,
    // where p is a 3-dimensional position, and t is a scalar.
    inline BoundVec3 pointAtParameter(const double t) const noexcept {
        return this->origin_ + (this->direction_ * t);
    }

    // Returns the time of intersection at a Point p.
    // The calculation: t = [p.a() - ray.origin().a()] / ray.unitDirection().a()
    //                  where a is a non-zero unit direction of the ray.
    // To reduce a vector multiplication to a single multiplication for the given direction,
    // we can do the following:
    // Since Point p = ray.origin() + ray.direction() * (v +/- discriminant),
    // We can simply provide the difference or addition of v and the discriminant.
    inline double timeOfIntersectionAt(double discriminant_v) const noexcept {
      return (NZD_params_.unit_direction * discriminant_v) * NZD_params_.inv_direction;
    }

    // Similar to above implementation, but uses a given vector p.
    inline double timeOfIntersectionAt(const Vec3 &p) const noexcept {
      return (p[NZD_params_.index] - NZD_params_.origin) * NZD_params_.inv_direction;
    }

    inline BoundVec3 origin() const noexcept { return this->origin_; }

    inline FreeVec3 direction() const noexcept { return this->direction_; }

    inline FreeVec3 invDirection() const noexcept { return this->inverse_direction_; }

    inline UnitVec3 unitDirection() const noexcept { return this->unit_direction_; }

private:
    // The origin of the ray.
    const BoundVec3 origin_;

    // The direction of the ray.
    const FreeVec3 direction_;

    // The normalized direction of the ray.
    const UnitVec3 unit_direction_;

    // The inverse direction of the ray.
    const FreeVec3 inverse_direction_;

    // Non-zero direction parameters.
    NonZeroDirectionParameters NZD_params_;
};

#endif //SPHERICAL_VOLUME_RENDERING_RAY_H
