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
    Vec3(const double x, const double y, const double z)
            : x_(x), y_(y), z_(z) {}

    inline double x() const { return this->x_; }

    inline double y() const { return this->y_; }

    inline double z() const { return this->z_; }

    inline double &x() { return this->x_; }

    inline double &y() { return this->y_; }

    inline double &z() { return this->z_; }

    inline double length() const {
        return std::sqrt(this->x() * this->x() + this->y() * this->y() + this->z() * this->z());
    }

    inline double squared_length() const {
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
struct FreeVec3 : Vec3 {

    explicit FreeVec3(const Vec3 &vec3) : Vec3(vec3.x(), vec3.y(), vec3.z()) {}

    explicit FreeVec3(double x, double y, double z) : Vec3(x, y, z) {}

    inline double dot(const Vec3 &other) const {
        return this->x() * other.x() + this->y() * other.y() + this->z() * other.z();
    }

    inline FreeVec3 &operator+=(const FreeVec3 &other) {
        this->x() += other.x();
        this->y() += other.y();
        this->z() += other.z();
        return *this;
    }

    inline FreeVec3 &operator-=(const FreeVec3 &other) {
        this->x() -= other.x();
        this->y() -= other.y();
        this->z() -= other.z();
        return *this;
    }

    inline FreeVec3 &operator*=(const double scalar) {
        this->x() *= scalar;
        this->y() *= scalar;
        this->z() *= scalar;
        return *this;
    }

    inline FreeVec3 &operator/=(const double scalar) {
        assert(scalar != 0.0);
        this->x() /= scalar;
        this->y() /= scalar;
        this->z() /= scalar;
        return *this;
    }
};

inline FreeVec3 operator+(const FreeVec3 &v) { return v; }

inline FreeVec3 operator-(const FreeVec3 &v) {
    return FreeVec3(-v.x(), -v.y(), -v.z());
}

inline FreeVec3 operator+(FreeVec3 v1, const FreeVec3 &v2) {
    return v1 += v2;
}

inline FreeVec3 operator-(FreeVec3 v1, const FreeVec3 &v2) {
    return v1 -= v2;
}

inline FreeVec3 operator*(FreeVec3 v, const double scalar) {
    return v *= scalar;
}

inline FreeVec3 operator/(FreeVec3 v, const double scalar) {
    return v /= scalar;
}

// A 3-dimensional bounded vector has a fixed start and end point. It represents a fixed point
// in space, relative to some frame of reference.
struct BoundVec3 : Vec3 {
    explicit BoundVec3(const Vec3 &vec3) : Vec3(vec3.x(), vec3.y(), vec3.z()) {}

    explicit BoundVec3(double x, double y, double z) : Vec3(x, y, z) {}

    inline double dot(const Vec3 &other) const {
        return this->x() * other.x() + this->y() * other.y() + this->z() * other.z();
    }

    inline BoundVec3 &operator+=(const FreeVec3 &other) {
        this->x() += other.x();
        this->y() += other.y();
        this->z() += other.z();
        return *this;
    }

    inline BoundVec3 &operator-=(const FreeVec3 &other) {
        return *this += (-other);
    }
};

inline FreeVec3 operator-(const BoundVec3 &v1, const BoundVec3 &v2) {
    return FreeVec3(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

inline BoundVec3 operator+(BoundVec3 v1, const FreeVec3 &v2) {
    return v1 += v2;
}

inline BoundVec3 operator-(BoundVec3 v1, const FreeVec3 &v2) {
    return v1 -= v2;
}

// Represents a 3-dimensional unit vector, an abstraction over free vectors that guarantees
// a length of 1. To prevent its length from changing, UnitVec3 does not allow
// for mutations.
struct UnitVec3 {
    explicit UnitVec3(double x, double y, double z)
            : UnitVec3(FreeVec3(x, y, z)) {}

    explicit UnitVec3(const Vec3 &vec3) : UnitVec3(FreeVec3(vec3)) {}

    explicit UnitVec3(const FreeVec3 &free_vec3) :
            inner_(free_vec3 / free_vec3.length()) {}

    inline double x() const { return this->to_free().x(); }

    inline double y() const { return this->to_free().y(); }

    inline double z() const { return this->to_free().z(); }

    inline const FreeVec3 &to_free() const { return inner_; }

private:
    const FreeVec3 inner_;
};

inline FreeVec3 operator*(const UnitVec3 &v, const double scalar) {
    return v.to_free() * scalar;
}

inline FreeVec3 operator/(const UnitVec3 &v, const double scalar) {
    assert(scalar != 0.0);
    return v.to_free() / scalar;
}

#endif //SPHERICAL_VOLUME_RENDERING_VEC3_H
