#include <algorithm>

#include "../spherical_volume_rendering_util.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// Testing for the spherical coordinate voxel traversal algorithm.
// Utilizes the Google Test suite.
//
// To build tests:
//   Install CMake version 3.7 or higher (https://cmake.org/). Then, run the
//   following command: >    cd cpp/tests && mkdir build && cd build && cmake ..
//   && make
// To run tests:
//   >    cd .. && ./bin/test_svr
//
// For information on Google Test, see:
// https://github.com/google/googletest/blob/master/googletest/README.md For
// examples of Google Test, see:
// https://github.com/google/googletest/tree/master/googletest/samples

namespace {
constexpr double TAU = 2 * M_PI;
constexpr svr::SphereBound MIN_BOUND = {
    .radial = 0.0, .polar = 0.0, .azimuthal = 0.0};

// Determines equality amongst actual spherical voxels, and the expected
// spherical voxels.
void verifyEqualVoxels(const std::vector<svr::SphericalVoxel> &actual_voxels,
                       const std::vector<int> &expected_radial_voxels,
                       const std::vector<int> &expected_theta_voxels,
                       const std::vector<int> &expected_phi_voxels) {
  const std::size_t num_voxels = actual_voxels.size();
  std::vector<int> radial_voxels(num_voxels);
  std::vector<int> theta_voxels(num_voxels);
  std::vector<int> phi_voxels(num_voxels);
  std::transform(
      actual_voxels.cbegin(), actual_voxels.cend(), radial_voxels.begin(),
      [](const svr::SphericalVoxel &voxel) -> int { return voxel.radial; });
  std::transform(
      actual_voxels.cbegin(), actual_voxels.cend(), theta_voxels.begin(),
      [](const svr::SphericalVoxel &voxel) -> int { return voxel.polar; });
  std::transform(
      actual_voxels.cbegin(), actual_voxels.cend(), phi_voxels.begin(),
      [](const svr::SphericalVoxel &voxel) -> int { return voxel.azimuthal; });
  EXPECT_THAT(radial_voxels, testing::ContainerEq(expected_radial_voxels));
  EXPECT_THAT(theta_voxels, testing::ContainerEq(expected_theta_voxels));
  EXPECT_THAT(phi_voxels, testing::ContainerEq(expected_phi_voxels));
}

TEST(SphericalCoordinateTraversal, RayDoesNotEnterSphere) {
  const BoundVec3 sphere_center(15.0, 15.0, 15.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 8;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(3.0, 3.0, 3.0);
  const UnitVec3 ray_direction(-2.0, -1.3, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  EXPECT_EQ(actual_voxels.size(), 0);
}

TEST(SphericalCoordinateTraversal, RayDoesNotEnterSphereTangentialHit) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 8;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-10.0, -10.0, 0.0);
  const UnitVec3 ray_direction(0.0, 1.0, 0.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  EXPECT_EQ(actual_voxels.size(), 0);
}

TEST(SphericalCoordinateTraversal, RayBeginsWithinSphere) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-3.0, 4.0, 5.0);
  const UnitVec3 ray_direction(1.0, -1.0, -1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {2, 3, 4, 4, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 0, 3, 3, 3, 3, 3};
  const std::vector<int> expected_phi_voxels = {1, 1, 1, 0, 0, 3, 3, 3, 3};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(DISABLED_SphericalCoordinateTraversal, RayEndsWithinSphere) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(13.0, -15.0, 16.0);
  const UnitVec3 ray_direction(-1.5, 1.2, -1.5);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 2, 3};
  const std::vector<int> expected_theta_voxels = {3, 3, 2, 2};
  const std::vector<int> expected_phi_voxels = {0, 0, 1, 1};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(DISABLED_SphericalCoordinateTraversal, RayBeginsAndEndsWithinSphere) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-3.0, 4.0, 5.0);
  const UnitVec3 ray_direction(1.0, -1.0, -1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {2, 3, 4, 4, 4};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 0, 3};
  const std::vector<int> expected_phi_voxels = {1, 1, 1, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(DISABLED_SphericalCoordinateTraversal,
     RayBeginsAndEndsWithinSphereNotCenteredAtOrigin) {
  const BoundVec3 sphere_center(2.0, 3.0, 2.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-1.0, 7.0, 7.0);
  const UnitVec3 ray_direction(1.0, -1.0, -1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {2, 3, 4, 4, 4};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 0, 3};
  const std::vector<int> expected_phi_voxels = {1, 1, 1, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, SphereCenteredAtOrigin) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-13.0, -13.0, -13.0);
  const UnitVec3 ray_direction(1.0, 1.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, SphereNotCenteredAtOrigin) {
  const BoundVec3 sphere_center(2.0, 2.0, 2.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-11.0, -11.0, -11.0);
  const UnitVec3 ray_direction(1.0, 1.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RaySlightOffsetInXYPlane) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-13.0, -13.0, -13.0);
  const UnitVec3 ray_direction(1.0, 1.5, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 2, 3, 2, 2, 1};
  const std::vector<int> expected_theta_voxels = {2, 2, 1, 1, 1, 0, 0};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 2, 2, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayTravelsAlongXAxis) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 8;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-15.0, 0.0, 0.0);
  const UnitVec3 ray_direction(1.0, 0.0, 0.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {3, 3, 3, 3, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {1, 1, 1, 1, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayTravelsAlongYAxis) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 8;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(0.0, -15.0, 0.0);
  const UnitVec3 ray_direction(0.0, 1.0, 0.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {5, 5, 5, 5, 1, 1, 1, 1};
  const std::vector<int> expected_phi_voxels = {0, 0, 0, 0, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayTravelsAlongZAxis) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 8;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(0.0, 0.0, -15.0);
  const UnitVec3 ray_direction(0.0, 0.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {0, 0, 0, 0, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayParallelToXYPlane) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-15.0, -15.0, 0.0);
  const UnitVec3 ray_direction(1.0, 1.0, 0.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {1, 1, 1, 1, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayParallelToXZPlane) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-15.0, 0.0, -15.0);
  const UnitVec3 ray_direction(1.0, 0.0, 1.0);
  const Ray ray(ray_origin, ray_direction);
  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 1, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayParallelToYZPlane) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(0.0, -15.0, -15.0);
  const UnitVec3 ray_direction(0.0, 1.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayDirectionNegativeXPositiveYZ) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(13.0, -15.0, -15.0);
  const UnitVec3 ray_direction(-1.0, 1.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {3, 3, 3, 2, 2, 1, 1, 1, 1};
  const std::vector<int> expected_phi_voxels = {3, 3, 3, 2, 2, 1, 1, 1, 1};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayDirectionNegativeYPositiveXZ) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-13.0, 17.0, -15.0);
  const UnitVec3 ray_direction(1.0, -1.2, 1.3);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 3, 4,
                                                   4, 3, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 1, 1, 0, 0, 3, 3, 3};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 1, 1, 0, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayDirectionNegativeZPositiveXY) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-13.0, -12.0, 15.3);
  const UnitVec3 ray_direction(1.4, 2.0, -1.3);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 1, 2, 2, 1};
  const std::vector<int> expected_theta_voxels = {2, 1, 1, 0, 0};
  const std::vector<int> expected_phi_voxels = {1, 1, 1, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayDirectionNegativeXYZ) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(15.0, 12.0, 15.0);
  const UnitVec3 ray_direction(-1.4, -2.0, -1.3);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 1, 2, 1, 1};
  const std::vector<int> expected_theta_voxels = {0, 3, 3, 3, 2};
  const std::vector<int> expected_phi_voxels = {0, 0, 0, 0, 1};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, OddNumberAngularSections) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 9.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 3;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-15.0, -15.0, -15.0);
  const UnitVec3 ray_direction(1.0, 1.0, 1.3);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 2, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 1, 0, 0};
  const std::vector<int> expected_phi_voxels = {2, 2, 1, 1, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, OddNumberAzimuthalSections) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 3;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-15.0, -15.0, -15.0);
  const UnitVec3 ray_direction(1.0, 1.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {1, 1, 1, 1, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, LargeNumberOfRadialSections) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 40;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-15.0, -15.0, -15.0);
  const UnitVec3 ray_direction(1.0, 1.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {
      1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16,
      17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
      33, 34, 35, 36, 37, 38, 39, 40, 40, 39, 38, 37, 36, 35, 34, 33,
      32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17,
      16, 15, 14, 13, 12, 11, 10, 9,  8,  7,  6,  5,  4,  3,  2,  1};
  const std::vector<int> expected_theta_voxels = {
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, LargeNumberOfAngularSections) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 40;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-15.0, -15.0, -15.0);
  const UnitVec3 ray_direction(1.0, 1.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {24, 24, 24, 24, 4, 4, 4, 4};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, LargeNumberOfAzimuthalSections) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 40;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-15.0, -15.0, -15.0);
  const UnitVec3 ray_direction(1.0, 1.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4, 4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {2, 2, 2, 2, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {24, 24, 24, 24, 4, 4, 4, 4};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(DISABLED_SphericalCoordinateTraversal,
     RayBeginsInOutermostRadiusAndEndsWithinSphere) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-4.0, -4.0, -6.0);
  const UnitVec3 ray_direction(1.3, 1.0, 1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 3, 4, 4};
  const std::vector<int> expected_theta_voxels = {2, 2, 2, 3, 3, 0};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 3, 3, 3};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayBeginsAtSphereOrigin) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(0.0, 0.0, 0.0);
  const UnitVec3 ray_direction(-1.5, 1.2, -1.5);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {4, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 1};
  const std::vector<int> expected_phi_voxels = {2, 2, 2, 2};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayBeginsPastSphereOriginOne) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-3.0, 2.4, -3.0);
  const UnitVec3 ray_direction(-1.5, 1.2, -1.5);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {3, 2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1, 1};
  const std::vector<int> expected_phi_voxels = {2, 2, 2};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayBeginsPastSphereOriginTwo) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-4.5, 3.6, -4.5);
  const UnitVec3 ray_direction(-1.5, 1.2, -1.5);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1};
  const std::vector<int> expected_phi_voxels = {2, 2};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayBeginsPastSphereOriginThree) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-6.0, 4.8, -6.0);
  const UnitVec3 ray_direction(-1.5, 1.2, -1.5);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1};
  const std::vector<int> expected_theta_voxels = {1};
  const std::vector<int> expected_phi_voxels = {2};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, RayBeginsPastSphereOriginFour) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-7.5, 6.0, -7.5);
  const UnitVec3 ray_direction(-1.5, 1.2, -1.5);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {};
  const std::vector<int> expected_theta_voxels = {};
  const std::vector<int> expected_phi_voxels = {};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, TangentialHitWithInnerRadialVoxelOne) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-5.0, 0.0, 10.0);
  const UnitVec3 ray_direction(0.0, 0.0, -1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 1};
  const std::vector<int> expected_phi_voxels = {1, 1, 2, 2};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, TangentialHitWithInnerRadialVoxelTwo) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-2.5, 0.0, 10.0);
  const UnitVec3 ray_direction(0.0, 0.0, -1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 1, 1, 1};
  const std::vector<int> expected_phi_voxels = {1, 1, 1, 2, 2, 2};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal,
     TangentialHitNoDoubleIntersectionWithSameVoxel) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 1;
  const std::size_t num_azimuthal_sections = 1;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-2.5, 0.0, 10.0);
  const UnitVec3 ray_direction(0.0, 0.0, -1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {0, 0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {0, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, NearlyTangentialHit) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = TAU};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const BoundVec3 ray_origin(-5.01, 0.0, 10.0);
  const UnitVec3 ray_direction(0.0, 0.0, -1.0);
  const Ray ray(ray_origin, ray_direction);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 2, 1};
  const std::vector<int> expected_theta_voxels = {1, 1, 1, 1};
  const std::vector<int> expected_phi_voxels = {1, 1, 2, 2};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(SphericalCoordinateTraversal, UpperHemisphereHit) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 8;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = M_PI};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(
      Ray(BoundVec3(-11.0, 2.0, 1.0), UnitVec3(1.0, 0.0, 0.0)), grid, t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 3, 4, 4,
                                                   4, 4, 3, 3, 2, 1};
  const std::vector<int> expected_theta_voxels = {3, 3, 3, 2, 2, 2,
                                                  1, 1, 1, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {3, 3, 3, 3, 3, 2,
                                                1, 0, 0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);

  const std::vector<BoundVec3> ray_origins = {
      BoundVec3(-5.0, -5.0, 5.0), BoundVec3(-1.0, -1.0, 10.0),
      BoundVec3(0.0, 0.0, 15.0), BoundVec3(-3.0, -3.0, 1.0),
      BoundVec3(-1.0, -5.0, 20.0)};
  for (const auto &ray_origin : ray_origins) {
    const UnitVec3 ray_direction(0.0, 0.0, -1.0);
    const auto v =
        walkSphericalVolume(Ray(ray_origin, ray_direction), grid, t_end);
    EXPECT_NE(v.size(), 0);
  }
}

TEST(SphericalCoordinateTraversal, UpperHemisphereMiss) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 8;
  const std::size_t num_azimuthal_sections = 4;
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = TAU, .azimuthal = M_PI};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);

  const double t_end = 1.0;
  const BoundVec3 ray_origin = BoundVec3(-5.0, -5.0, -5.0);
  const UnitVec3 ray_direction(1.0, 0.0, 0.0);
  const auto actual_voxels =
      walkSphericalVolume(Ray(ray_origin, ray_direction), grid, t_end);
  EXPECT_EQ(actual_voxels.size(), 0);
}

TEST(SphericalCoordinateTraversal, AvoidRaySteppingToRadialVoxelZero) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10e3;
  const std::size_t num_radial_sections = 128;
  const std::size_t num_polar_sections = 128;
  const std::size_t num_azimuthal_sections = 128;
  const svr::SphereBound min_bound = {
      .radial = 0.0, .polar = 0.0, .azimuthal = 0.0};
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = 2 * M_PI, .azimuthal = 2 * M_PI};
  const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);

  const double t_end = 1.0;
  const Ray ray(BoundVec3(-984.375, 250.0, -10001.0), UnitVec3(0.0, 0.0, 1.0));
  const auto actual_voxels = walkSphericalVolume(ray, grid, t_end);
  const std::size_t last_idx = actual_voxels.size() - 1;
  EXPECT_NE(actual_voxels[last_idx].radial, 0);
}

TEST(SphericalCoordinateTraversal, VerifyManyRaysEntranceAndExit) {
  // Given an orthographic ray projection with sufficient time, all rays should
  // enter and exit the sphere.
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10e4;
  const std::size_t num_radial_sections = 32;
  const std::size_t num_polar_sections = 32;
  const std::size_t num_azimuthal_sections = 32;
  const svr::SphereBound min_bound = {
      .radial = 0.0, .polar = 0.0, .azimuthal = 0.0};
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = 2 * M_PI, .azimuthal = 2 * M_PI};
  const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);

  const double t_end = 1.0;
  {  // Z-axis
    const UnitVec3 ray_direction(0.0, 0.0, 1.0);
    double ray_origin_x = -1000.0;
    double ray_origin_y = -1000.0;
    const double ray_origin_z = -(sphere_max_radius + 1.0);

    const double ray_origin_plane_movement = 2000.0 / 30;
    for (std::size_t i = 0; i < 30; ++i) {
      for (std::size_t j = 0; j < 30; ++j) {
        const BoundVec3 ray_origin(ray_origin_x, ray_origin_y, ray_origin_z);
        const auto actual_voxels =
            walkSphericalVolume(Ray(ray_origin, ray_direction), grid, t_end);
        EXPECT_NE(actual_voxels.size(), 0);
        EXPECT_EQ(actual_voxels[0].radial, 1);
        const std::size_t last = actual_voxels.size() - 1;
        EXPECT_EQ(actual_voxels[last].radial, 1);
        ray_origin_y =
            (j == 30 - 1) ? -1000.0 : ray_origin_y + ray_origin_plane_movement;
      }
      ray_origin_x += ray_origin_plane_movement;
    }
  }
  {  // Y-axis
    const UnitVec3 ray_direction(0.0, 1.0, 0.0);
    double ray_origin_x = -1000.0;
    double ray_origin_z = -1000.0;
    const double ray_origin_y = -(sphere_max_radius + 1.0);

    const double ray_origin_plane_movement = 2000.0 / 30;
    for (std::size_t i = 0; i < 30; ++i) {
      for (std::size_t j = 0; j < 30; ++j) {
        const BoundVec3 ray_origin(ray_origin_x, ray_origin_y, ray_origin_z);
        const auto actual_voxels =
            walkSphericalVolume(Ray(ray_origin, ray_direction), grid, t_end);
        EXPECT_NE(actual_voxels.size(), 0);
        EXPECT_EQ(actual_voxels[0].radial, 1);
        const std::size_t last = actual_voxels.size() - 1;
        EXPECT_EQ(actual_voxels[last].radial, 1);
        ray_origin_z =
            (j == 30 - 1) ? -1000.0 : ray_origin_z + ray_origin_plane_movement;
      }
      ray_origin_x += ray_origin_plane_movement;
    }
  }
  {  // X-axis
    const UnitVec3 ray_direction(1.0, 0.0, 0.0);
    double ray_origin_y = -1000.0;
    double ray_origin_z = -1000.0;
    const double ray_origin_x = -(sphere_max_radius + 1.0);

    const double ray_origin_plane_movement = 2000.0 / 30;
    for (std::size_t i = 0; i < 30; ++i) {
      for (std::size_t j = 0; j < 30; ++j) {
        const BoundVec3 ray_origin(ray_origin_x, ray_origin_y, ray_origin_z);
        const auto actual_voxels =
            walkSphericalVolume(Ray(ray_origin, ray_direction), grid, t_end);
        EXPECT_NE(actual_voxels.size(), 0);
        EXPECT_EQ(actual_voxels[0].radial, 1);
        const std::size_t last = actual_voxels.size() - 1;
        EXPECT_EQ(actual_voxels[last].radial, 1);
        ray_origin_y =
            (j == 30 - 1) ? -1000.0 : ray_origin_y + ray_origin_plane_movement;
      }
      ray_origin_z += ray_origin_plane_movement;
    }
  }
}

TEST(DISABLED_SphericalCoordinateTraversal, FirstQuadrantHit) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 8;
  const svr::SphereBound max_bound = {.radial = sphere_max_radius,
                                      .polar = M_PI / 2.0,
                                      .azimuthal = M_PI / 2.0};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(
      Ray(BoundVec3(13.0, 13.0, 13.0), UnitVec3(-1.0, -1.0, -1.0)), grid,
      t_end);
  const std::vector<int> expected_radial_voxels = {1, 2, 3, 4};
  const std::vector<int> expected_theta_voxels = {0, 0, 0, 0};
  const std::vector<int> expected_phi_voxels = {0, 0, 0, 0};
  verifyEqualVoxels(actual_voxels, expected_radial_voxels,
                    expected_theta_voxels, expected_phi_voxels);
}

TEST(DISABLED_SphericalCoordinateTraversal, FirstQuadrantMiss) {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10.0;
  const std::size_t num_radial_sections = 4;
  const std::size_t num_polar_sections = 4;
  const std::size_t num_azimuthal_sections = 8;
  const svr::SphereBound max_bound = {.radial = sphere_max_radius,
                                      .polar = M_PI / 2.0,
                                      .azimuthal = M_PI / 2.0};
  const svr::SphericalVoxelGrid grid(MIN_BOUND, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);

  const double t_end = 1.0;
  const auto actual_voxels = walkSphericalVolume(
      Ray(BoundVec3(13.0, -13.0, 13.0), UnitVec3(-1.0, 1.0, -1.0)), grid,
      t_end);
  EXPECT_EQ(actual_voxels.size(), 0);
}

}  // namespace