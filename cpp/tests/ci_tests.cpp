#include <algorithm>
#include <random>

#include "../spherical_volume_rendering_util.h"
#include "gtest/gtest.h"

// A set of tests to be run by continuous integration. These verify basic
// traversal properties such as ordering and bounds.
// todo: This was written hastily and likely can be refactored to reduce
//       time complexity and code duplication.
namespace {
void printRayData(const Ray& ray) {
  printf("\nRay origin: {%f, %f, %f}", ray.origin().x(), ray.origin().y(),
         ray.origin().z());
  printf("\nRay direction: {%f, %f, %f}", ray.direction().x(),
         ray.direction().y(), ray.direction().z());
}

void printVoxelInformation(const svr::SphericalVoxel v,
                           const std::string& info = "") {
  if (info.empty()) {
    printf("{%d, %d, %d} ", v.radial, v.polar, v.azimuthal);
  } else {
    printf("\nAbout: %s\n Voxel: {%d, %d, %d}\n", info.data(), v.radial,
           v.polar, v.azimuthal);
  }
}

template <class It>
void printVoxelGroupingInformation(const std::vector<svr::SphericalVoxel>& v,
                                   It it, const std::string& info) {
  printf("\n%s\n", info.data());
  if (it - 1 == v.cbegin()) {
    printf("\n[Voxel is the first]");
  } else {
    printVoxelInformation(*(it - 1));
  }
  printVoxelInformation(*it);
  if (it + 1 == v.cend()) {
    printf("\n[Voxel is the last]");
  } else {
    printVoxelInformation(*(it + 1));
  }
}

// Verifies the entrance and radial voxel is 1 for all rays, and each radial
// voxel is within bounds (0, number_of_radial_voxels]. Also verifies
// the number of voxels for each traversal is greater than two. Lastly,
// verifies each radial voxel's transition order. If its solely a radial hit,
// the next radial voxel should be current+1 or current-1. If it is not, then
// the next radial voxel is current+1, current-1, or current+0.
bool CheckRadialVoxelsForOrthographicProjection(
    const Ray& ray, const std::vector<svr::SphericalVoxel>& actual_voxels,
    int number_of_radial_voxels) {
  const auto it = std::adjacent_find(
      actual_voxels.cbegin(), actual_voxels.cend(),
      [](const svr::SphericalVoxel& v1, const svr::SphericalVoxel& v2) {
        const bool radial_hit_only =
            v1.polar == v2.polar && v1.azimuthal == v2.azimuthal;
        if (radial_hit_only) {
          return !(v1.radial - 1 == v2.radial || v1.radial + 1 == v2.radial);
        }
        const bool within_one =
            (v1.radial == v2.radial || v1.radial - 1 == v2.radial ||
             v1.radial + 1 == v2.radial);
        return !within_one;
      });
  if (it != actual_voxels.cend()) {
    printVoxelGroupingInformation(
        actual_voxels, it,
        "The current radial voxel is not within +- 1 of the next voxel.");
    printRayData(ray);
    return false;
  }
  if (actual_voxels.empty()) {
    printf("\nNo intersection with sphere at all.");
    printRayData(ray);
    EXPECT_FALSE(actual_voxels.empty());
    return false;
  }
  const std::size_t last = actual_voxels.size() - 1;
  if (actual_voxels.size() < 2) {
    printf("\nRay traversed less than two voxels.");
    printRayData(ray);
    EXPECT_GE(actual_voxels.size(), 2);
    return false;
  }
  if (actual_voxels[0].radial != 1 || actual_voxels[last].radial != 1) {
    printf("\nDid not complete entire traversal.");
    const auto first_voxel = actual_voxels[0];
    const auto last_voxel = actual_voxels[last];
    printRayData(ray);
    printVoxelInformation(first_voxel, "Entrance Voxel.");
    printVoxelInformation(last_voxel, "Exit Voxel");
    EXPECT_TRUE(actual_voxels[0].radial == 1);
    EXPECT_TRUE(actual_voxels[last].radial == 1);
    return false;
  }
  const auto it2 = std::find_if_not(
      actual_voxels.cbegin(), actual_voxels.cend(), [&](svr::SphericalVoxel i) {
        return 0 < i.radial <= number_of_radial_voxels;
      });
  if (it2 != actual_voxels.cend()) {
    printRayData(ray);
    printVoxelGroupingInformation(
        actual_voxels, it2,
        "There exists a radial voxel i such that"
        "0 < i <= number_of_radial_voxels does not hold.");
    return false;
  }
  return true;
}

// It should hold true in orthographic projects that each angular voxel should
// be within +- 1 of the last angular voxel except for at most in one case. This
// case occurs when traversing the line x = 0.
bool checkOrthographicProjectionAngularVoxelOrdering(
    const Ray& ray, const std::vector<svr::SphericalVoxel>& v) {
  const auto it_polar = std::adjacent_find(
      v.cbegin(), v.cend(),
      [](const svr::SphericalVoxel& v1, const svr::SphericalVoxel& v2) {
        const bool polar_within_one =
            (v1.polar == v2.polar || v1.polar - 1 == v2.polar ||
             v1.polar + 1 == v2.polar);
        return !polar_within_one;
      });
  if (it_polar != v.cend()) {
    const auto it2_polar = std::adjacent_find(
        it_polar + 1, v.cend(),
        [](const svr::SphericalVoxel& v1, const svr::SphericalVoxel& v2) {
          const bool polar_within_one =
              (v1.polar == v2.polar || v1.polar - 1 == v2.polar ||
               v1.polar + 1 == v2.polar);
          return !polar_within_one;
        });
    if (it2_polar != v.cend()) {
      printRayData(ray);
      printVoxelGroupingInformation(
          v, it2_polar,
          "A polar voxel makes two jumps greater than +-1 "
          "voxel. This should only occur once per ray when the ray passes "
          "the line X = 0.");
      printVoxelGroupingInformation(v, it_polar, "Previous Jump:");
      return false;
    }
  }
  const auto it_azimuthal = std::adjacent_find(
      v.cbegin(), v.cend(),
      [](const svr::SphericalVoxel& v1, const svr::SphericalVoxel& v2) {
        const bool azimuthal_within_one =
            (v1.azimuthal == v2.azimuthal || v1.azimuthal - 1 == v2.azimuthal ||
             v1.azimuthal + 1 == v2.azimuthal);
        return !azimuthal_within_one;
      });
  if (it_azimuthal != v.cend()) {
    const auto it2_azimuthal = std::adjacent_find(
        it_azimuthal + 1, v.cend(),
        [](const svr::SphericalVoxel& v1, const svr::SphericalVoxel& v2) {
          const bool azimuthal_within_one = (v1.azimuthal == v2.azimuthal ||
                                             v1.azimuthal - 1 == v2.azimuthal ||
                                             v1.azimuthal + 1 == v2.azimuthal);
          return !azimuthal_within_one;
        });
    if (it2_azimuthal != v.cend()) {
      printRayData(ray);
      printVoxelGroupingInformation(
          v, it2_azimuthal,
          "An azimuthal voxel makes two jumps greater than +-1 "
          "voxel. This should only occur once per ray when the ray passes "
          "the line X = 0.");
      printVoxelGroupingInformation(v, it_azimuthal, "Previous Jump:");
      return false;
    }
  }
  return true;
}

// Verifies all polar and azimuthal voxels are within bounds
// 0 <= X <= number_of_sections.
bool CheckAngularVoxels(const std::vector<svr::SphericalVoxel>& v,
                        const Ray& ray, int number_of_angular_sections) {
  const auto it =
      std::find_if_not(v.cbegin(), v.cend(), [&](svr::SphericalVoxel i) {
        return 0 <= i.azimuthal && i.azimuthal < number_of_angular_sections &&
               0 <= i.polar && i.polar < number_of_angular_sections;
      });
  if (it != v.cend()) {
    printRayData(ray);
    printVoxelInformation(
        *it,
        "\n There exists an angular voxel i such that"
        "0 <= i <= number_of_angular_sections does not hold.");
    return false;
  }
  return checkOrthographicProjectionAngularVoxelOrdering(ray, v);
}

// Sends X^2 rays through a Y^3 spherical voxel grid orthographically. All rays
// are perpendicular to the Z plane.
void inline orthographicTraverseXSquaredRaysinYCubedVoxels(
    const std::size_t X, const std::size_t Y) noexcept {
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10e4;
  const std::size_t num_radial_sections = Y;
  const std::size_t num_polar_sections = Y;
  const std::size_t num_azimuthal_sections = Y;
  const svr::SphereBound min_bound = {
      .radial = 0.0, .polar = 0.0, .azimuthal = 0.0};
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = 2 * M_PI, .azimuthal = 2 * M_PI};
  const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const double t_begin = 0.0;
  const double t_end = sphere_max_radius * 3;

  const FreeVec3 ray_direction(0.0, 0.0, 1.0);
  double ray_origin_x = -1000.0;
  double ray_origin_y = -1000.0;
  const double ray_origin_z = -(sphere_max_radius + 1.0);

  const double ray_origin_plane_movement = 2000.0 / X;
  for (std::size_t i = 0; i < X; ++i) {
    for (std::size_t j = 0; j < X; ++j) {
      const BoundVec3 ray_origin(ray_origin_x, ray_origin_y, ray_origin_z);
      const Ray ray(ray_origin, ray_direction);
      const auto actual_voxels = walkSphericalVolume(ray, grid, t_begin, t_end);
      if (!CheckRadialVoxelsForOrthographicProjection(ray, actual_voxels, Y)) {
        const bool radial_voxels_check_passed = false;
        EXPECT_TRUE(radial_voxels_check_passed);
        return;
      }
      if (!CheckAngularVoxels(actual_voxels, ray, Y)) {
        const bool angular_voxels_check_passed = false;
        EXPECT_TRUE(angular_voxels_check_passed);
        return;
      }
      ray_origin_y =
          (j == X - 1) ? -1000.0 : ray_origin_y + ray_origin_plane_movement;
    }
    ray_origin_x += ray_origin_plane_movement;
  }
}

// Similar to orthographicTraverseXSquaredRaysinYCubedVoxels, but uses a
// seeded random direction within bounds [1.0, 3.0), and random origin
// within bounds [0.0, 10,000.0). All rays are perpendicular to the y-axis.
// The number of sections for each voxel type is bounded by [16, Y];
void inline randomRayPlacementTraverseXSquaredRaysInYBoundedCubedVoxels(
    const std::size_t X, const std::size_t Y) noexcept {
  std::default_random_engine rd;
  std::mt19937 mt(rd());
  EXPECT_GT(Y, 24);
  std::uniform_int_distribution<int> num_sections(16, Y);
  const BoundVec3 sphere_center(0.0, 0.0, 0.0);
  const double sphere_max_radius = 10e6;
  const std::size_t num_radial_sections = num_sections(mt);
  const std::size_t num_polar_sections = num_sections(mt);
  const std::size_t num_azimuthal_sections = num_sections(mt);
  const svr::SphereBound min_bound = {
      .radial = 0.0, .polar = 0.0, .azimuthal = 0.0};
  const svr::SphereBound max_bound = {
      .radial = sphere_max_radius, .polar = 2 * M_PI, .azimuthal = 2 * M_PI};
  const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                     num_polar_sections, num_azimuthal_sections,
                                     sphere_center);
  const double t_begin = 0.0;
  const double t_end = sphere_max_radius * 100;

  const double ray_origin_y = -(sphere_max_radius + 1.0);
  for (int i = 0; i < X * X; ++i) {
    std::uniform_real_distribution<double> dist1(0.0, 10000.0);
    const double ray_origin_x = dist1(mt);
    const double ray_origin_z = dist1(mt);
    const BoundVec3 ray_origin(ray_origin_x, ray_origin_y, ray_origin_z);
    std::uniform_real_distribution<double> dist2(1.0, 3.0);
    const double ray_x_dir = dist2(mt);
    const double ray_y_dir = dist2(mt);
    const double ray_z_dir = dist2(mt);
    const Ray ray(ray_origin, FreeVec3(ray_x_dir, ray_y_dir, ray_z_dir));
    const auto actual_voxels = walkSphericalVolume(ray, grid, t_begin, t_end);
    if (!CheckRadialVoxelsForOrthographicProjection(ray, actual_voxels, Y)) {
      const bool radial_voxels_check_passed = false;
      EXPECT_TRUE(radial_voxels_check_passed);
      return;
    }
    if (!CheckAngularVoxels(actual_voxels, ray, Y)) {
      const bool angular_voxels_check_passed = false;
      EXPECT_TRUE(angular_voxels_check_passed);
      return;
    }
  }
}

struct TestParameters {
  std::size_t ray_squared_count;
  std::size_t voxel_cubed_count;
};

const std::vector<TestParameters> random_test_parameters = {
    {.ray_squared_count = 32, .voxel_cubed_count = 32},
    {.ray_squared_count = 64, .voxel_cubed_count = 32},
    {.ray_squared_count = 64, .voxel_cubed_count = 64},
    {.ray_squared_count = 128, .voxel_cubed_count = 64},
    {.ray_squared_count = 64, .voxel_cubed_count = 128},
    {.ray_squared_count = 128, .voxel_cubed_count = 128},
};

const std::vector<TestParameters> orthographic_test_parameters = {
    {.ray_squared_count = 64, .voxel_cubed_count = 64},
    {.ray_squared_count = 128, .voxel_cubed_count = 64},
    {.ray_squared_count = 256, .voxel_cubed_count = 64},
    {.ray_squared_count = 64, .voxel_cubed_count = 128},
    {.ray_squared_count = 128, .voxel_cubed_count = 128},
    {.ray_squared_count = 64, .voxel_cubed_count = 512},
    {.ray_squared_count = 64, .voxel_cubed_count = 1024},
    {.ray_squared_count = 512, .voxel_cubed_count = 32},
    {.ray_squared_count = 1024, .voxel_cubed_count = 32},
};

TEST(ContinuousIntegration, RandomRayPlacement) {
  for (const auto param : random_test_parameters) {
    printf("   [ RUN      ] %lu^2 Rays in [16, %lu^3] Voxels\n",
           param.ray_squared_count, param.voxel_cubed_count);
    randomRayPlacementTraverseXSquaredRaysInYBoundedCubedVoxels(
        param.ray_squared_count, param.voxel_cubed_count);
    printf("   [       OK ] %lu^2 Rays in [16, %lu^3] Voxels\n",
           param.ray_squared_count, param.voxel_cubed_count);
  }
}

TEST(ContinuousIntegration, OrthographicProjection) {
  for (const auto param : orthographic_test_parameters) {
    printf("   [ RUN      ] %lu^2 Rays in %lu^3 Voxels\n",
           param.ray_squared_count, param.voxel_cubed_count);
    orthographicTraverseXSquaredRaysinYCubedVoxels(param.ray_squared_count,
                                                   param.voxel_cubed_count);
    printf("   [       OK ] %lu^2 Rays in %lu^3 Voxels\n",
           param.ray_squared_count, param.voxel_cubed_count);
  }
}

}  // namespace
