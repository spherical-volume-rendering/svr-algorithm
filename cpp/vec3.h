#ifndef SPHERICAL_VOLUME_RENDERING_VEC3_H
#define SPHERICAL_VOLUME_RENDERING_VEC3_H

#include <algorithm>
#include <numeric>
#include <cmath>
#include <array>

// The indices for Vec3. For example, Vec3[0] returns the x-direction.
enum DirectionIndex {
    X_DIRECTION = 0,
    Y_DIRECTION = 1,
    Z_DIRECTION = 2
};

// Represents a Euclidean vector in 3-dimensional space.
// Assumes vectors take the form of:
//      [x]
//      [y]
//      [z]
struct Vec3 {
public:
    constexpr inline Vec3(const double x, const double y, const double z) : e_{x,y,z} {}

    constexpr inline double x() const noexcept { return this->e_[0]; }

    constexpr inline double y() const noexcept { return this->e_[1]; }

    constexpr inline double z() const noexcept { return this->e_[2]; }

    inline double &x() noexcept { return this->e_[0]; }

    inline double &y() noexcept { return this->e_[1]; }

    inline double &z() noexcept { return this->e_[2]; }

    inline double length() const noexcept {
        return std::sqrt(this->e_[0] * this->e_[0]  + this->e_[1]  * this->e_[1]  + this->e_[2]  * this->e_[2]);
    }

    constexpr inline double squared_length() const noexcept {
        return e_[0] * e_[0] + e_[1] * e_[1] + e_[2] * e_[2];
    }

    inline double operator[](const std::size_t index) const noexcept {
        return e_[index];
  }

private:
    std::array<double, 3> e_;
};

// A 3-dimensional free vector, which has no initial point. It has two main criteria:
// (1) direction, and (2) magnitude.
// TODO(cgyurgyik): While code duplication will occur, it may be more efficient to avoid inheritance.
//                  This will eliminate VTABLE lookups, though this may be (a) negligible or (b) non-existent. See:
//                  Langr, J. "Modern C++ Programming with Test-Driven Development: Code Better, Sleep Better" [5.10]
struct FreeVec3 : Vec3 {

    constexpr inline explicit FreeVec3(const Vec3 &vec3) : Vec3(vec3.x(), vec3.y(), vec3.z()) {}

    constexpr inline explicit FreeVec3(double x, double y, double z) : Vec3(x, y, z) {}

    constexpr inline double dot(const Vec3 &other) const noexcept {
        return this->x() * other.x() + this->y() * other.y() + this->z() * other.z();
    }

    inline FreeVec3 &operator+=(const FreeVec3 &other) noexcept {
        this->x() += other.x();
        this->y() += other.y();
        this->z() += other.z();
        return *this;
    }

    inline FreeVec3 &operator-=(const FreeVec3 &other) noexcept {
        this->x() -= other.x();
        this->y() -= other.y();
        this->z() -= other.z();
        return *this;
    }

    inline FreeVec3 &operator*=(const double scalar) noexcept {
        this->x() *= scalar;
        this->y() *= scalar;
        this->z() *= scalar;
        return *this;
    }

    inline FreeVec3 &operator/=(const double scalar) noexcept {
        this->x() /= scalar;
        this->y() /= scalar;
        this->z() /= scalar;
        return *this;
    }
};

inline FreeVec3 operator+(const FreeVec3 &v) noexcept { return v; }

inline FreeVec3 operator-(const FreeVec3 &v) noexcept {
    return FreeVec3(-v.x(), -v.y(), -v.z());
}

inline FreeVec3 operator+(FreeVec3 v1, const FreeVec3 &v2) noexcept {
    return v1 += v2;
}

inline FreeVec3 operator-(FreeVec3 v1, const FreeVec3 &v2) noexcept {
    return v1 -= v2;
}

inline FreeVec3 operator*(FreeVec3 v, const double scalar) noexcept {
    return v *= scalar;
}

inline FreeVec3 operator/(FreeVec3 v, const double scalar) noexcept {
    return v /= scalar;
}

// A 3-dimensional bounded vector has a fixed start and end point. It represents a fixed point
// in space, relative to some frame of reference.
struct BoundVec3 : Vec3 {
    constexpr inline explicit BoundVec3(const Vec3 &vec3) : Vec3(vec3.x(), vec3.y(), vec3.z()) {}

    constexpr inline explicit BoundVec3(double x, double y, double z) : Vec3(x, y, z) {}

    constexpr inline double dot(const Vec3 &other) const noexcept {
        return this->x() * other.x() + this->y() * other.y() + this->z() * other.z();
    }

    inline BoundVec3 &operator+=(const FreeVec3 &other) noexcept {
        this->x() += other.x();
        this->y() += other.y();
        this->z() += other.z();
        return *this;
    }

    inline BoundVec3 &operator-=(const FreeVec3 &other) noexcept {
        return *this += (-other);
    }
};

inline FreeVec3 operator-(const BoundVec3 &v1, const BoundVec3 &v2) noexcept {
    return FreeVec3(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

inline BoundVec3 operator+(BoundVec3 v1, const FreeVec3 &v2) noexcept {
    return v1 += v2;
}

inline BoundVec3 operator-(BoundVec3 v1, const FreeVec3 &v2) noexcept {
    return v1 -= v2;
}

// Represents a 3-dimensional unit vector, an abstraction over free vectors that guarantees
// a length of 1. To prevent its length from changing, UnitVec3 does not allow
// for mutations.
struct UnitVec3 {
    inline explicit UnitVec3(double x, double y, double z)
            : UnitVec3(FreeVec3(x, y, z)) {}

    inline explicit UnitVec3(const Vec3 &vec3) : UnitVec3(FreeVec3(vec3)) {}

    inline explicit UnitVec3(const FreeVec3 &free_vec3) :
            inner_(free_vec3 / free_vec3.length()) {}

    inline double x() const noexcept { return this->to_free().x(); }

    inline double y() const noexcept { return this->to_free().y(); }

    inline double z() const noexcept { return this->to_free().z(); }

    inline const FreeVec3 &to_free() const noexcept { return inner_; }

    inline double operator[](const std::size_t index) const noexcept {
        return this->to_free()[index];
    }

private:
    const FreeVec3 inner_;
};

inline FreeVec3 operator*(const UnitVec3 &v, const double scalar) noexcept {
    return v.to_free() * scalar;
}

inline FreeVec3 operator/(const UnitVec3 &v, const double scalar) noexcept {
    return v.to_free() / scalar;
}

#endif //SPHERICAL_VOLUME_RENDERING_VEC3_H
