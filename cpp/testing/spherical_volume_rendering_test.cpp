#include "gtest/gtest.h"
#include "../spherical_volume_rendering_util.h"

// Utilizes the Google Test suite.
// To run:
//       1. Clone the google test suite to spherical-volume-rendering/cpp/testing
//          - It can be found at: https://github.com/google/googletest.git
//          - Currently uses the default folder name (googletest) in the CMake file, so there is no need to change this.
//       2. Load (or re-load) CMakeLists.txt
//
// For information on Google Test, see: https://github.com/google/googletest/blob/master/googletest/README.md
// For examples of Google Test, see: https://github.com/google/googletest/tree/master/googletest/samples
namespace {
    // The tolerance provided for the floating point error carried during the
    // spherical coordinate traversal algorithm.
    constexpr double tolerance = 10e-20;

    TEST(SphericalCoordinateTraversal, RayDoesNotEnterSphere) {
        const Vec3 min_bound(0.0, 0.0, 0.0);
        const Vec3 max_bound(30.0, 30.0, 30.0);
        const Vec3 sphere_center(15.0, 15.0, 15.0);
        const double sphere_max_radius = 10.0;
        const size_t num_radial_sections = 4;
        const size_t num_angular_sections = 8;
        const size_t num_azimuthal_sections = 4;
        const SphericalVoxelGrid grid(BoundVec3(min_bound), BoundVec3(max_bound), num_radial_sections,
                                      num_angular_sections,
                                      num_azimuthal_sections, BoundVec3(sphere_center), sphere_max_radius);
        const Vec3 ray_origin(3.0, 3.0, 3.0);
        const UnitVec3 ray_direction(-2.0, -1.3, 1.0);
        const Ray ray(BoundVec3(ray_origin), ray_direction);

        const double t_begin = 0.0;
        const double t_end = 15.0;
        const auto voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end, tolerance);
        EXPECT_EQ(voxels.size(), 0);
    }
}