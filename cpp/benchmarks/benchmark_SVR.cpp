#include <benchmark/benchmark.h>

#include "../spherical_volume_rendering_util.h"

#define DEBUG 0

// Benchmarking for the spherical coordinate voxel traversal algorithm.
// Utilises the Google Benchmark library.
//
// To build the benchmarks:
//  Install CMake version 3.7 or higher (https://cmake.org/). Then, run the
//  following command: >    cd cpp/benchmarks && mkdir build && cd build &&
//  cmake .. && make
//
//  To run the benchmarks:
//  >    cd .. && ./bin/benchmark_SVR
//
// For more information on Google Benchmark, see:
// https://github.com/google/benchmark For examples of Google Benchmark, see:
// https://github.com/google/benchmark#usage

namespace {

// Sends X^2 rays through a Y^3 voxel sphere with maximum radius 10e4.
// The set up is the following:
// This traversal is orthographic in nature, and all rays will intersect the
// sphere. In the X plane: Ray origin moves incrementally from [-1,000.0,
// 1,000.0]. In the Y plane: Ray origin moves incrementally from [-1,000.0,
// 1,000.0]. In the Z plane: Ray origin begins at -(10e4 + 1.0). It does not
// move in the Z plane. From this, one can infer that the ray moves
// incrementally in the XY plane from
// (-1,000.0, -1,000.0) -> (1,000.0, 1,000.0) while remaining outside the sphere
// in the Z plane. Since the maximum sphere radius is 10e4, this ensures all
// rays will intersect.
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
      const auto actual_voxels = walkSphericalVolume(
          Ray(ray_origin, ray_direction), grid, t_begin, t_end);
#if DEBUG
      if (actual_voxels.size() == 0) {
        printf("\n No intersection at all.");
      } else {
        const std::size_t last = actual_voxels.size() - 1;
        if (actual_voxels[0].radial != 1 || actual_voxels[last].radial != 1 ||
            actual_voxels.size() < 2) {
          printf("\nDid not complete entire traversal.");
          const auto first_voxel = actual_voxels[0];
          const auto last_voxel = actual_voxels[last];
          printf("\nRay origin: {%f, %f, %f}", ray_origin_x, ray_origin_y,
                 ray_origin_z);
          printf("\nEntrance Voxel: {%d, %d, %d} ... Exit Voxel: {%d, %d, %d}",
                 first_voxel.radial, first_voxel.polar, first_voxel.azimuthal,
                 last_voxel.radial, last_voxel.polar, last_voxel.azimuthal);
          for (const auto v : actual_voxels)
            printf("{%d, %d, %d} , ", v.radial, v.polar, v.azimuthal);
        }
      }
#endif
      ray_origin_y =
          (j == X - 1) ? -1000.0 : ray_origin_y + ray_origin_plane_movement;
    }
    ray_origin_x += ray_origin_plane_movement;
  }
}

static void Orthographic_128SquaredRays_64CubedVoxels(benchmark::State &state) {
  for (auto _ : state) {
    orthographicTraverseXSquaredRaysinYCubedVoxels(128, 64);
  }
}

static void Orthographic_256SquaredRays_64CubedVoxels(benchmark::State &state) {
  for (auto _ : state) {
    orthographicTraverseXSquaredRaysinYCubedVoxels(256, 64);
  }
}

static void Orthographic_512SquaredRays_64CubedVoxels(benchmark::State &state) {
  for (auto _ : state) {
    orthographicTraverseXSquaredRaysinYCubedVoxels(512, 64);
  }
}

static void Orthographic_128SquaredRays_128CubedVoxels(
    benchmark::State &state) {
  for (auto _ : state) {
    orthographicTraverseXSquaredRaysinYCubedVoxels(128, 128);
  }
}

static void Orthographic_256SquaredRays_128CubedVoxels(
    benchmark::State &state) {
  for (auto _ : state) {
    orthographicTraverseXSquaredRaysinYCubedVoxels(256, 128);
  }
}

static void Orthographic_512SquaredRays_128CubedVoxels(
    benchmark::State &state) {
  for (auto _ : state) {
    orthographicTraverseXSquaredRaysinYCubedVoxels(512, 128);
  }
}

BENCHMARK(Orthographic_128SquaredRays_64CubedVoxels)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(1);
BENCHMARK(Orthographic_256SquaredRays_64CubedVoxels)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(1);
BENCHMARK(Orthographic_512SquaredRays_64CubedVoxels)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(1);
BENCHMARK(Orthographic_128SquaredRays_128CubedVoxels)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(1);
BENCHMARK(Orthographic_256SquaredRays_128CubedVoxels)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(1);
BENCHMARK(Orthographic_512SquaredRays_128CubedVoxels)
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(1);

}  // namespace

BENCHMARK_MAIN();
