#ifndef SPHERICAL_VOLUME_RENDERING_RAY_H
#define SPHERICAL_VOLUME_RENDERING_RAY_H

#include "vec3.h"

// Encapsulates the functionality of a ray.
// This consists of two components, the origin of the ray,
// and the direction of the ray.
struct Ray final {
     Ray(const BoundVec3& origin, const FreeVec3& direction)
            : origin_(origin), direction_(direction), unit_direction_(direction),
            inverse_direction_(FreeVec3(1.0 / direction.x(), 1.0 / direction.y(), 1.0 / direction.z())),
              x_dir_is_non_zero_(std::abs(direction.x()) > 0.0), y_dir_is_non_zero_(std::abs(direction.y()) > 0.0),
              z_dir_is_non_zero_(std::abs(direction.z()) > 0.0) {}

    // Represents the function p(t) = origin + t * direction,
    // where p is a 3-dimensional position, and t is a scalar.
    inline BoundVec3 pointAtParameter(const double t) const {
        return this->origin_ + (this->direction_ * t);
    }

    // Returns the time of intersection at a Point p.
    // The calculation: t = [p.a() - ray.origin().a()] / ray.unitDirection().a()
    //                  where a is a non-zero unit direction of the ray.
    // To reduce a vector multiplication to a single multiplication for the given direction,
    // we can do the following:
    // Since Point p = ray.origin() + ray.direction() * (v +/- discriminant),
    // We can simply provide the difference or addition of v and the discriminant.
    inline double timeOfIntersectionAt(double discriminant_v) const {
         if (xDirectionIsNonZero()) {
             const double p_x = origin().x() + unitDirection().x() * discriminant_v;
             return (p_x - origin().x()) * invDirection().x();
         }
         if (yDirectionIsNonZero()) {
             const double p_y = origin().y() + unitDirection().y() * discriminant_v;
             return (p_y - origin().y()) * invDirection().y();
         }
        const double p_z = origin().z() + unitDirection().z() * discriminant_v;
        return (p_z - origin().z()) * invDirection().z();
     }

    inline BoundVec3 origin() const { return this->origin_; }
    inline FreeVec3 direction() const { return this->direction_; }
    inline FreeVec3 invDirection() const { return this->inverse_direction_; }
    inline UnitVec3 unitDirection() const { return this->unit_direction_; }
    inline bool xDirectionIsNonZero() const { return x_dir_is_non_zero_; }
    inline bool yDirectionIsNonZero() const { return y_dir_is_non_zero_; }
    inline bool zDirectionIsNonZero() const { return z_dir_is_non_zero_; }

private:
    // The origin of the ray.
    const BoundVec3 origin_;
    // The direction of the ray.
    const FreeVec3 direction_;
    // The normalized direction of the ray.
    const UnitVec3 unit_direction_;
    // The inverse direction of the ray.
    const FreeVec3 inverse_direction_;
    // Determines whether the given direction is non-zero to avoid division by 0.
    const bool x_dir_is_non_zero_, y_dir_is_non_zero_, z_dir_is_non_zero_;
};

#endif //SPHERICAL_VOLUME_RENDERING_RAY_H
