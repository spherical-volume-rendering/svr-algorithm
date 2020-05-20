#include <algorithm>

#include "../spherical_volume_rendering_util.h"
#include "gtest/gtest.h"

// A set of tests to be run by continuous integration. These test general
// functionality of sending many rays through a spherical volume.
// todo: This was written hastily and likely can be refactored to reduce
//       time complexity and code duplication.
namespace {

// Verifies the entrance and radial voxel is 1 for all rays, and each radial
// voxel is within bounds (0, number_of_radial_voxels]. Also verifies
// the number of voxels for each traversal is greater than two. Lastly,
// verifies each radial voxel is within +-1 of the last radial voxel.
bool CheckRadialVoxelsForOrthographicProjection(
    const Ray& ray, const std::vector<svr::SphericalVoxel>& actual_voxels,
    int number_of_radial_voxels) {
  const auto it = std::adjacent_find(
      actual_voxels.cbegin(), actual_voxels.cend(),
      [](const svr::SphericalVoxel& v1, const svr::SphericalVoxel& v2) {
        const bool within_one =
            (v1.radial == v2.radial || v1.radial - 1 == v2.radial ||
             v1.radial + 1 == v2.radial);
        return !within_one;
      });
  if (it != actual_voxels.cend()) {
    const auto v = *it;
    printf("\nThe current radial voxel is not within +- 1 of the next voxel.");
    printf("Current Voxel: {%d, %d, %d}", v.radial, v.polar, v.azimuthal);
    if (it - 1 == actual_voxels.cbegin()) {
      printf("\n Previous Voxel: None. This is the first voxel");
    } else {
      const auto vb = *(it - 1);
      printf("\nPrevious Voxel: {%d, %d, %d}", vb.radial, vb.polar,
             vb.azimuthal);
    }
    if (it + 1 == actual_voxels.cend()) {
      printf("\n Next Voxel: None. This is the last voxel");
    } else {
      const auto vn = *(it + 1);
      printf("\nNext Voxel: {%d, %d, %d}", vn.radial, vn.polar, vn.azimuthal);
    }
    printf("\nRay origin: {%f, %f, %f}", ray.origin().x(), ray.origin().y(),
           ray.origin().z());
    return false;
  }
  if (actual_voxels.size() == 0) {
    printf("\n No intersection at all.");
    printf("\nRay origin: {%f, %f, %f}", ray.origin().x(), ray.origin().y(),
           ray.origin().z());
    EXPECT_FALSE(actual_voxels.size() == 0);
    return false;
  }
  const std::size_t last = actual_voxels.size() - 1;
  if (actual_voxels.size() < 2) {
    printf("\nRay traversed less than two voxels.");
    printf("\nRay origin: {%f, %f, %f}", ray.origin().x(), ray.origin().y(),
           ray.origin().z());
    EXPECT_GE(actual_voxels.size(), 2);
    return false;
  }
  if (actual_voxels[0].radial != 1 || actual_voxels[last].radial != 1) {
    printf("\nDid not complete entire traversal.");
    const auto first_voxel = actual_voxels[0];
    const auto last_voxel = actual_voxels[last];
    printf("\nRay origin: {%f, %f, %f}", ray.origin().x(), ray.origin().y(),
           ray.origin().z());
    printf("\nEntrance Voxel: {%d, %d, %d} ... Exit Voxel: {%d, %d, %d}",
           first_voxel.radial, first_voxel.polar, first_voxel.azimuthal,
           last_voxel.radial, last_voxel.polar, last_voxel.azimuthal);
    EXPECT_TRUE(actual_voxels[0].radial == 1);
    EXPECT_TRUE(actual_voxels[last].radial == 1);
    return false;
  }
  const auto it2 = std::find_if_not(
      actual_voxels.cbegin(), actual_voxels.cend(), [&](svr::SphericalVoxel i) {
        return 0 < i.radial <= number_of_radial_voxels;
      });
  if (it2 != actual_voxels.cend()) {
    const auto vxl = *it2;
    printf(
        "\n There exists a radial voxel i such that"
        "0 < i <= number_of_radial_voxels does not hold.");
    printf("\nRay origin: {%f, %f, %f}", ray.origin().x(), ray.origin().y(),
           ray.origin().z());
    printf("\nVoxel: {%d, %d, %d}", vxl.radial, vxl.polar, vxl.azimuthal);
    return false;
  }
  return true;
}

// It should hold that each angular voxel should be within +- 1 of the last
// angular voxel except for at most in one case. This case occurs
// when traversing the line x = 0.
bool checkAngularVoxelOrdering(const Ray& ray,
                                const std::vector<svr::SphericalVoxel>& v,
                                int number_of_angular_voxels) {
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
        it_polar, v.cend(),
        [](const svr::SphericalVoxel& v1, const svr::SphericalVoxel& v2) {
          const bool polar_within_one =
              (v1.polar == v2.polar || v1.polar - 1 == v2.polar ||
               v1.polar + 1 == v2.polar);
          return !polar_within_one;
        });
    if (it2_polar != v.cend()) {
      const auto vxl = *it2_polar;
      printf(
          "\n A polar voxel makes two jumps greater than +-1 voxel. This"
          "should only occur once per ray when the ray passes the line X = 0.");
      printf("\nRay origin: {%f, %f, %f}", ray.origin().x(), ray.origin().y(),
             ray.origin().z());
      printf("\nVoxel: {%d, %d, %d}", vxl.radial, vxl.polar, vxl.azimuthal);
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
        it_polar, v.cend(),
        [](const svr::SphericalVoxel& v1, const svr::SphericalVoxel& v2) {
          const bool azimuthal_within_one = (v1.azimuthal == v2.azimuthal ||
                                             v1.azimuthal - 1 == v2.azimuthal ||
                                             v1.azimuthal + 1 == v2.azimuthal);
          return !azimuthal_within_one;
        });
    if (it2_azimuthal != v.cend()) {
      const auto vxl = *it2_azimuthal;
      printf(
          "\n An azimuthal voxel makes two jumps greater than +-1 voxel. This"
          "should only occur once per ray when the ray passes the line X = 0.");
      printf("\nRay origin: {%f, %f, %f}", ray.origin().x(), ray.origin().y(),
             ray.origin().z());
      printf("\nVoxel: {%d, %d, %d}", vxl.radial, vxl.polar, vxl.azimuthal);
      return false;
    }
  }
  return true;
}

// Verifies all polar and azimuthal voxels are within bounds
// 0 <= X <= number_of_voxels.
bool CheckAngularVoxelsForOrthographicProjection(
    const std::vector<svr::SphericalVoxel>& v, const Ray& ray,
    int number_of_angular_voxels) {
  const auto it =
      std::find_if_not(v.cbegin(), v.cend(), [&](svr::SphericalVoxel i) {
        return 0 <= i.azimuthal && i.azimuthal < number_of_angular_voxels &&
               0 <= i.polar && i.polar < number_of_angular_voxels;
      });
  if (it != v.cend()) {
    const auto vxl = *it;
    printf(
        "\n There exists an angular voxel i such that"
        "0 <= i <= number_of_angular_voxels does not hold.");
    printf("\nRay origin: {%f, %f, %f}", ray.origin().x(), ray.origin().y(),
           ray.origin().z());
    printf("\nVoxel: {%d, %d, %d}", vxl.radial, vxl.polar, vxl.azimuthal);
    return false;
  }
  return checkAngularVoxelOrdering(ray, v, number_of_angular_voxels);
}

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
      const auto actual_voxels = walkSphericalVolume(
          Ray(ray_origin, ray_direction), grid, t_begin, t_end);
      if (!CheckRadialVoxelsForOrthographicProjection(ray, actual_voxels, Y)) {
        return;
      }
      if (!CheckAngularVoxelsForOrthographicProjection(actual_voxels, ray, Y)) {
        return;
      }
      ray_origin_y =
          (j == X - 1) ? -1000.0 : ray_origin_y + ray_origin_plane_movement;
    }
    ray_origin_x += ray_origin_plane_movement;
  }
}

TEST(ContinuousIntegration, 64SquaredRaysIn64CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(64, 64);
}

TEST(ContinuousIntegration, 128SquaredRaysIn64CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(128, 64);
}

TEST(ContinuousIntegration, 256SquaredRaysIn64CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(256, 64);
}

TEST(ContinuousIntegration, 512SquaredRaysIn64CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(512, 64);
}

TEST(ContinuousIntegration, 1024SquaredRaysIn64CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(1024, 64);
}

TEST(ContinuousIntegration, 64SquaredRaysIn128CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(64, 128);
}

TEST(ContinuousIntegration, 128SquaredRaysIn128CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(128, 128);
}

TEST(ContinuousIntegration, 64SquaredRaysIn256CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(64, 256);
}

TEST(ContinuousIntegration, 64SquaredRaysIn512CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(64, 512);
}

TEST(ContinuousIntegration, 64SquaredRaysIn1024CubedVoxels) {
  orthographicTraverseXSquaredRaysinYCubedVoxels(64, 1024);
}

}  // namespace
