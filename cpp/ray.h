#ifndef SPHERICAL_VOLUME_RENDERING_RAY_H
#define SPHERICAL_VOLUME_RENDERING_RAY_H

#include "vec3.h"

// Determines the non-zero direction for a given ray direction.
constexpr DirectionIndex getNonZeroDirection(
    const UnitVec3 &direction) noexcept {
  return direction.x() != 0.0
             ? X_DIRECTION
             : direction.y() != 0.0 ? Y_DIRECTION : Z_DIRECTION;
}

// Calculates the inverse of each component of a unit direction.
constexpr FreeVec3 inverseDirection(const UnitVec3 &direction) {
  return FreeVec3(1.0 / direction.x(), 1.0 / direction.y(),
                  1.0 / direction.z());
}

// Encapsulates the functionality of a ray. This consists of two components, the
// origin of the ray, and the unit direction of the ray. To avoid checking for a
// non-zero direction upon each function call, these parameters are initialized
// upon construction.
struct Ray final {
  inline Ray(const BoundVec3 &origin, const UnitVec3 &direction) noexcept
      : origin_(origin),
        direction_(direction),
        inverse_direction_(inverseDirection(direction)),
        NZD_index_(getNonZeroDirection(direction)) {}

  // Represents the function p(t) = origin + t * direction,
  // where p is a 3-dimensional position, and t is a scalar.
  inline BoundVec3 pointAtParameter(const double t) const noexcept {
    return this->origin_ + this->direction_ * t;
  }

  // Returns the time of intersection at a Point p.
  // The calculation: t = [p.a() - ray.origin().a()] / ray.direction().a()
  //                  where a is a non-zero direction of the ray.
  // To reduce a vector multiplication to a single multiplication for the given
  // direction, we can do the following: Since Point p = ray.origin() +
  // ray.direction() * (v +/- discriminant), We can simply provide the
  // difference or addition of v and the discriminant.
  inline double timeOfIntersectionAt(double discriminant_v) const noexcept {
    return this->direction_[NZD_index_] * discriminant_v *
           this->inverse_direction_[NZD_index_];
  }

  // Similar to above implementation, but uses a given vector p.
  inline double timeOfIntersectionAt(const Vec3 &p) const noexcept {
    return (p[NZD_index_] - this->origin_[NZD_index_]) *
           this->inverse_direction_[NZD_index_];
  }

  inline const BoundVec3 &origin() const noexcept { return this->origin_; }

  inline const UnitVec3 &direction() const noexcept { return this->direction_; }

  inline const FreeVec3 &invDirection() const noexcept {
    return this->inverse_direction_;
  }

  inline DirectionIndex NonZeroDirectionIndex() const noexcept {
    return this->NZD_index_;
  }

 private:
  // The origin of the ray.
  const BoundVec3 origin_;

  // The direction of the ray.
  const UnitVec3 direction_;

  // The inverse direction of the ray.
  const FreeVec3 inverse_direction_;

  // Index of a non-zero direction.
  const enum DirectionIndex NZD_index_;
};

// Pre-calculated information for the generalized angular hit function, which
// generalizes azimuthal and polar hits. Since the ray segment is dependent
// solely on time, this is unnecessary to calculate twice for each plane hit
// function. Here, ray_segment is the difference between P2 and P1.
struct RaySegment {
 public:
  inline RaySegment(double max_t, const Ray &ray) noexcept
      : P2_(ray.pointAtParameter(max_t)), NZDI_(ray.NonZeroDirectionIndex()) {}

  // Updates the point P1 with the new time traversal time t. Similarly, updates
  // the segment denoted by P2 - P1.
  inline void updateAtTime(double t, const Ray &ray) noexcept {
    P1_ = ray.pointAtParameter(t);
    ray_segment_ = P2_ - P1_;
  }

  // Calculates the updated ray segment intersection point given an intersect
  // parameter. More information on the use case can be found at:
  // http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
  inline double intersectionTimeAt(double intersect_parameter,
                                   const Ray &ray) const noexcept {
    return (P1_[NZDI_] + ray_segment_[NZDI_] * intersect_parameter -
            ray.origin()[NZDI_]) *
           ray.invDirection()[NZDI_];
  }

  inline const BoundVec3 &P1() const noexcept { return P1_; }

  inline const BoundVec3 &P2() const noexcept { return P2_; }

  inline const FreeVec3 &vector() const noexcept { return ray_segment_; }

 private:
  // The end point of the ray segment.
  const BoundVec3 P2_;

  // The non-zero direction index of the ray.
  const DirectionIndex NZDI_;

  // The begin point of the ray segment.
  BoundVec3 P1_;

  // The free vector represented by P2 - P1.
  FreeVec3 ray_segment_;
};

#endif  // SPHERICAL_VOLUME_RENDERING_RAY_H
