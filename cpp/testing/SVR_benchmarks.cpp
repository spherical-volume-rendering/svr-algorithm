#include <benchmark/benchmark.h>
#include "../spherical_volume_rendering_util.h"

// Benchmarking for the Spherical Volume Rendering algorithm.
// Utilises the Google Benchmark library found at:
// https://github.com/google/benchmark
//
// To use Google Benchmark, see: https://github.com/google/benchmark#installation

namespace {
    // Sends X^2 rays through a Y^3 voxel sphere with maximum radius 10^4.
    // The set up is the following:
    // This traversal is orthographic in nature, and all rays will intersect the sphere.
    // In the X plane: Ray origin moves incrementally from [-1000.0, 1000.0].
    // In the Y plane: Ray origin moves incrementally from [-1000.0, 1000.0].
    // In the Z plane: Ray origin begins at -(10^4 + 1000.0). It does not move in the Z plane.
    // From this, one can infer that the ray moves incrementally in the XY plane from
    // (-1000.0, -1000.0) -> (1000.0, 1000.0) while remaining outside the sphere in the Z plane.
    // Since the maximum sphere radius is 10^4, this ensures all rays will intersect.
    void inline orthographicTraverseXSquaredRaysinYCubedVoxels(const std::size_t X, const std::size_t Y) noexcept {
        const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
        const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
        const BoundVec3 sphere_center(0.0, 0.0, 0.0);
        const double sphere_max_radius = 10.0 * 1000.0;
        const std::size_t num_radial_sections = Y;
        const std::size_t num_angular_sections = Y;
        const std::size_t num_azimuthal_sections = Y;
        const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections, num_angular_sections,
                                           num_azimuthal_sections, sphere_center, sphere_max_radius);
        const double t_begin = 0.0;
        const double t_end = sphere_max_radius * 2.0;

        double ray_origin_x = -1000.0;
        double ray_origin_y = -1000.0;
        const double ray_origin_z = -(sphere_max_radius + 1000.0);

        const double ray_origin_plane_movement = 2000.0 / X;
        for (std::size_t i = 0; i < X; ++i) {
            for (std::size_t j = 0; j < X; ++j) {
                const BoundVec3 ray_origin(ray_origin_x, ray_origin_y, ray_origin_z);
                const FreeVec3  ray_direction(0.0, 0.0, 1.0);
                const auto actual_voxels = sphericalCoordinateVoxelTraversal(Ray(ray_origin, ray_direction),
                                                                             grid, t_begin, t_end);
                ray_origin_y = (j == X - 1) ? -1000.0 : ray_origin_y + ray_origin_plane_movement;
            }
            ray_origin_x += ray_origin_plane_movement;
        }
    }

    static void TraversalOne(benchmark::State &state) {
        for (auto _ : state) {
            // Traversal from bottom left to upper right.
            const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center(0.0, 0.0, 0.0);
            const double sphere_max_radius = 10000.0;
            const std::size_t num_radial_sections = 10000;
            const std::size_t num_angular_sections = 10000;
            const std::size_t num_azimuthal_sections = 10000;
            const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                               num_angular_sections,
                                               num_azimuthal_sections, sphere_center, sphere_max_radius);
            const BoundVec3 ray_origin(-13000.0, -13000.0, -13000.0);
            const FreeVec3 ray_direction(1.0, 1.0, 1.0);
            const Ray ray(ray_origin, ray_direction);
            const double t_begin = 0.0;
            const double t_end = 100000.0;
            const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
        }
    }

    static void TraversalTwo(benchmark::State &state) {
        for (auto _ : state) {
            // Traversal from upper right to bottom left.
            const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center(0.0, 0.0, 0.0);
            const double sphere_max_radius = 10000.0;
            const std::size_t num_radial_sections = 10000;
            const std::size_t num_angular_sections = 10000;
            const std::size_t num_azimuthal_sections = 10000;
            const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                               num_angular_sections,
                                               num_azimuthal_sections, sphere_center, sphere_max_radius);
            const BoundVec3 ray_origin(13000.0, 13000.0, 13000.0);
            const FreeVec3 ray_direction(-1.0, -1.0, -1.0);
            const Ray ray(ray_origin, ray_direction);
            const double t_begin = 0.0;
            const double t_end = 100000.0;
            const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
        }
    }

    static void TraversalParallelX(benchmark::State &state) {
        for (auto _ : state) {
            const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center(0.0, 0.0, 0.0);
            const double sphere_max_radius = 10000.0;
            const std::size_t num_radial_sections = 10000;
            const std::size_t num_angular_sections = 10000;
            const std::size_t num_azimuthal_sections = 10000;
            const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                               num_angular_sections,
                                               num_azimuthal_sections, sphere_center, sphere_max_radius);
            const BoundVec3 ray_origin(-13000.0, 10.0, 10.0);
            const FreeVec3 ray_direction(1.0, 0.0, 0.0);
            const Ray ray(ray_origin, ray_direction);
            const double t_begin = 0.0;
            const double t_end = 100000.0;
            const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
        }
    }

    static void TraversalParallelY(benchmark::State &state) {
        for (auto _ : state) {
            const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center(0.0, 0.0, 0.0);
            const double sphere_max_radius = 10000.0;
            const std::size_t num_radial_sections = 10000;
            const std::size_t num_angular_sections = 10000;
            const std::size_t num_azimuthal_sections = 10000;
            const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                               num_angular_sections,
                                               num_azimuthal_sections, sphere_center, sphere_max_radius);
            const BoundVec3 ray_origin(10.0, -13000.0, 0.0);
            const FreeVec3 ray_direction(0.0, 1.0, 0.0);
            const Ray ray(ray_origin, ray_direction);
            const double t_begin = 0.0;
            const double t_end = 100000.0;
            const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
        }
    }

    static void TraversalParallelZ(benchmark::State &state) {
        for (auto _ : state) {
            const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center(0.0, 0.0, 0.0);
            const double sphere_max_radius = 10000.0;
            const std::size_t num_radial_sections = 10000;
            const std::size_t num_angular_sections = 10000;
            const std::size_t num_azimuthal_sections = 10000;
            const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                               num_angular_sections,
                                               num_azimuthal_sections, sphere_center, sphere_max_radius);
            const BoundVec3 ray_origin(10.0, 0.0, -13000.0);
            const FreeVec3 ray_direction(0.0, 0.0, 1.0);
            const Ray ray(ray_origin, ray_direction);
            const double t_begin = 0.0;
            const double t_end = 100000.0;
            const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
        }
    }

    static void MultipleRayNoIntersection(benchmark::State &state) {
        for (auto _ : state) {
            const int number_of_runs = 100000;
            const BoundVec3 min_bound(0.0, 0.0, 0.0);
            const BoundVec3 max_bound(30.0, 30.0, 30.0);
            const BoundVec3 sphere_center(15.0, 15.0, 15.0);
            const double sphere_max_radius = 10.0;
            const std::size_t num_radial_sections = 4;
            const std::size_t num_angular_sections = 8;
            const std::size_t num_azimuthal_sections = 4;
            const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                               num_angular_sections,
                                               num_azimuthal_sections, sphere_center, sphere_max_radius);
            const BoundVec3 ray_origin(3.0, 3.0, 3.0);
            const FreeVec3 ray_direction(-2.0, -1.3, 1.0);
            const Ray ray(ray_origin, ray_direction);
            const double t_begin = 0.0;
            const double t_end = 15.0;
            for (int i = 0; i < number_of_runs; ++i) {
                sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
            }
        }
    }

    static void Orthographic_128SquaredRays_64CubedVoxels(benchmark::State &state) {
        for (auto _ : state) {
            orthographicTraverseXSquaredRaysinYCubedVoxels(128, 64);
        }
    }

    static void Orthographic_256SquaredRays_128CubedVoxels(benchmark::State &state) {
        for (auto _ : state) {
            orthographicTraverseXSquaredRaysinYCubedVoxels(256, 128);
        }
    }

    BENCHMARK(TraversalOne)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalTwo)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalParallelX)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalParallelY)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalParallelZ)->Unit(benchmark::kMillisecond);
    BENCHMARK(MultipleRayNoIntersection)->Unit(benchmark::kMillisecond);
    BENCHMARK(Orthographic_128SquaredRays_64CubedVoxels)->Unit(benchmark::kMillisecond);
    BENCHMARK(Orthographic_256SquaredRays_128CubedVoxels)->Unit(benchmark::kMillisecond)->Repetitions(10);

} // namespace

BENCHMARK_MAIN();
