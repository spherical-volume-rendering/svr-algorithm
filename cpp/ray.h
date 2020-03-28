#ifndef SPHERICAL_VOLUME_RENDERING_RAY_H
#define SPHERICAL_VOLUME_RENDERING_RAY_H

#include "vec3.h"

// Encapsulates the functionality of a ray.
// This consists of two components, the origin of the ray,
// and the direction of the ray.
struct Ray final {
     Ray(const BoundVec3& origin, const UnitVec3& direction)
            : origin_(origin), direction_(direction),
            inverse_direction_(UnitVec3(1.0 / direction.x(), 1.0 / direction.y(), 1.0 / direction.z())) {}

    // Represents the function p(t) = origin + t * direction,
    // where p is a 3-dimensional position, and t is a scalar.
    inline BoundVec3 point_at_parameter(const double t) const {
        return this->origin_ + (this->direction_ * t);
    }

    inline BoundVec3 origin() const { return this->origin_; }
    inline UnitVec3 direction() const { return this->direction_; }
    inline UnitVec3 inv_direction() const { return this->inverse_direction_;}
private:
    // The origin of the ray.
    const BoundVec3 origin_;
    // The normalized direction of the ray.
    const UnitVec3 direction_;
    // The normalized inverse direction of the ray.
    const UnitVec3 inverse_direction_;
};

#endif //SPHERICAL_VOLUME_RENDERING_RAY_H
