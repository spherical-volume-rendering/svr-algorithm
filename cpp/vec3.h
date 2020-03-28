#ifndef SPHERICAL_VOLUME_RENDERING_VEC3_H
#define SPHERICAL_VOLUME_RENDERING_VEC3_H

#include <algorithm>
#include <numeric>
#include <cmath>

// Represents a Euclidean vector in 3-dimensional space.
// Assumes vectors take the form of:
//      [x]
//      [y]
//      [z]
struct Vec3 {
public:
    constexpr explicit Vec3(const double x, const double y, const double z)
            : x_{x}, y_{y}, z_{z} {}

    [[nodiscard]] inline constexpr double x() const { return this->x_; }

    [[nodiscard]] inline constexpr double y() const { return this->y_; }

    [[nodiscard]] inline constexpr double z() const { return this->z_; }

    [[nodiscard]] inline constexpr double &x() { return this->x_; }

    [[nodiscard]] inline constexpr double &y() { return this->y_; }

    [[nodiscard]] inline constexpr double &z() { return this->z_; }

    [[nodiscard]] inline double length() const {
        return std::sqrt(this->x() * this->x() + this->y() * this->y() + this->z() * this->z());
    }

    [[nodiscard]] inline double squared_length() const {
        return x() * x() + y() * y() + z() * z();
    }

private:
    // Represents the x-dimension value of the vector.
    double x_;
    // Represents the y-dimension value of the vector.
    double y_;
    // Represents the z-dimension value of the vector.
    double z_;
};

// A 3-dimensional free vector, which has no initial point. It has two main criteria:
// (1) direction, and (2) magnitude.
// TODO(cgyurgyik): While code duplication will occur, it may be more efficient to avoid inheritance.
//                  This will eliminate VTABLE lookups, though this may be (a) negligible or (b) non-existent. See:
//                  Langr, J. "Modern C++ Programming with Test-Driven Development: Code Better, Sleep Better" [5.10]
struct FreeVec3 final : Vec3 {
    using Vec3::Vec3;

    constexpr explicit FreeVec3(const Vec3 &vec3) : Vec3::Vec3{vec3} {}

    [[nodiscard]] inline constexpr double dot(const Vec3 &other) const {
        return this->x() * other.x() + this->y() * other.y() + this->z() * other.z();
    }

    [[nodiscard]] inline constexpr FreeVec3 cross(const Vec3 &other) const {
        return FreeVec3{this->y() * other.z() - this->z() * other.y(),
                        this->z() * other.x() - this->x() * other.z(),
                        this->x() * other.y() - this->y() * other.x()};
    }

    [[nodiscard]] inline constexpr FreeVec3 &operator+=(const FreeVec3 &other) {
        this->x() += other.x();
        this->y() += other.y();
        this->z() += other.z();
        return *this;
    }

    [[nodiscard]] inline constexpr FreeVec3 &operator-=(const FreeVec3 &other) {
        this->x() -= other.x();
        this->y() -= other.y();
        this->z() -= other.z();
        return *this;
    }

    [[nodiscard]]  inline constexpr FreeVec3 &operator*=(const double scalar) {
        this->x() *= scalar;
        this->y() *= scalar;
        this->z() *= scalar;
        return *this;
    }

    [[nodiscard]] inline constexpr FreeVec3 &operator/=(const double scalar) {
        this->x() /= scalar;
        this->y() /= scalar;
        this->z() /= scalar;
        return *this;
    }
};

[[nodiscard]] inline constexpr FreeVec3 operator+(const FreeVec3 &v) { return v; }

[[nodiscard]] inline constexpr FreeVec3 operator-(const FreeVec3 &v) {
    return FreeVec3{-v.x(), -v.y(), -v.z()};
}

[[nodiscard]] inline constexpr FreeVec3 operator+(FreeVec3 v1, const FreeVec3 &v2) {
    return v1 += v2;
}

[[nodiscard]] inline constexpr FreeVec3 operator-(FreeVec3 v1, const FreeVec3 &v2) {
    return v1 -= v2;
}

[[nodiscard]] inline constexpr FreeVec3 operator*(FreeVec3 v, const double scalar) {
    return v *= scalar;
}

[[nodiscard]] inline constexpr FreeVec3 operator/(FreeVec3 v, const double scalar) {
    return v /= scalar;
}

// A 3-dimensional bounded vector has a fixed start and end point. It represents a fixed point
// in space, relative to some frame of reference.
struct BoundVec3 final : Vec3 {
    constexpr explicit BoundVec3(const Vec3 &vec3) : Vec3::Vec3{vec3} {}

    [[nodiscard]] inline constexpr double dot(const Vec3 &other) const {
        return this->x() * other.x() + this->y() * other.y() + this->z() * other.z();
    }

    [[nodiscard]] inline constexpr BoundVec3 &operator+=(const FreeVec3 &other) {
        this->x() += other.x();
        this->y() += other.y();
        this->z() += other.z();
        return *this;
    }

    [[nodiscard]] inline constexpr BoundVec3 &operator-=(const FreeVec3 &other) {
        return *this += (-other);
    }
};

[[nodiscard]] inline constexpr FreeVec3 operator-(const BoundVec3 &v1, const BoundVec3 &v2) {
    return FreeVec3{v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z()};
}

[[nodiscard]] inline constexpr BoundVec3 operator+(BoundVec3 v1, const FreeVec3 &v2) {
    return v1 += v2;
}

[[nodiscard]] inline constexpr BoundVec3 operator-(BoundVec3 v1, const FreeVec3 &v2) {
    return v1 -= v2;
}

// Represents a 3-dimensional unit vector, an abstraction over free vectors that guarantees
// a length of 1. To prevent its length from changing, UnitVec3 does not allow
// for mutations.
struct UnitVec3 final {
    constexpr UnitVec3(double x, double y, double z)
            : UnitVec3{FreeVec3{x, y, z}} {}

    constexpr explicit UnitVec3(const Vec3 &vec3) : UnitVec3{FreeVec3{vec3}} {}

    constexpr explicit UnitVec3(const FreeVec3 &free_vec3) :
            inner_{free_vec3 / free_vec3.length()} {}

    [[nodiscard]] inline constexpr double x() const { return this->to_free().x(); }

    [[nodiscard]] inline constexpr double y() const { return this->to_free().y(); }

    [[nodiscard]] inline constexpr double z() const { return this->to_free().z(); }

    [[nodiscard]] inline constexpr const FreeVec3 &to_free() const { return inner_; }

private:
    const FreeVec3 inner_;
};

[[nodiscard]] inline constexpr FreeVec3 operator*(const UnitVec3 &v, const double scalar) {
    return v.to_free() * scalar;
}

[[nodiscard]] inline constexpr FreeVec3 operator/(const UnitVec3 &v, const double scalar) {
    return v.to_free() / scalar;
}

#endif //SPHERICAL_VOLUME_RENDERING_VEC3_H
