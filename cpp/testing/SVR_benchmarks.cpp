#include <benchmark/benchmark.h>
#include "../spherical_volume_rendering_util.h"

// Benchmarking for the Spherical Volume Rendering algorithm.
// Utilises the Google Benchmark library found at:
// https://github.com/google/benchmark
//
// To use Google Benchmark, see: https://github.com/google/benchmark#installation

namespace {

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

    // 128^3 domain with 256^2 rays for a scratch paper benchmark.
    // In this case, we can imagine the ray's form a cube around the sphere,
    // and thus not all rays will intersect.
    static void OrthographicRayTracing(benchmark::State &state) {
        for (auto _ : state) {
            const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center(0.0, 0.0, 0.0);
            const double sphere_max_radius = 1000.0;
            const std::size_t num_radial_sections = 128;
            const std::size_t num_angular_sections = 128;
            const std::size_t num_azimuthal_sections = 128;
            const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                               num_angular_sections, num_azimuthal_sections,
                                               sphere_center, sphere_max_radius);
            const double t_begin = 0.0;
            const double t_end = 10000.0;

            double ray_origin_x = -1000.0;
            double ray_origin_y = -1000.0;
            const double ray_origin_z = -1100.0;
            const double ray_movement = 2000.0 / 256.0;
            for (int i = 0; i < 256; ++i) {
                for (int j = 0; j < 256; ++j) {
                    const BoundVec3 ray_origin(ray_origin_x, ray_origin_y, ray_origin_z);
                    const FreeVec3 ray_direction(0.0, 0.0, 1.0);
                    const Ray ray(ray_origin, ray_direction);
                    const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
                    ray_origin_y = (j != 255) ? ray_origin_y + ray_movement : -1000.0;
                }
                ray_origin_x += ray_movement;
            }
        }
    }

    // 128^3 domain with 256^2 rays for a scratch paper benchmark.
    // In this case, all rays intersect the sphere.
    static void OrthographicAllRaysIntersect(benchmark::State &state) {
        for (auto _ : state) {
            const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center(0.0, 0.0, 0.0);
            const double sphere_max_radius = 1000.0 * 10.0;
            const std::size_t num_radial_sections = 128;
            const std::size_t num_angular_sections = 128;
            const std::size_t num_azimuthal_sections = 128;
            const svr::SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                               num_angular_sections, num_azimuthal_sections,
                                               sphere_center, sphere_max_radius);
            const double t_begin = 0.0;
            const double t_end = 100000.0;

            double ray_origin_x = -1000.0;
            double ray_origin_y = -1000.0;
            const double ray_origin_z = -sphere_max_radius - 100.0;
            const double ray_movement = 2000.0 / 256.0;
            for (int i = 0; i < 256; ++i) {
                for (int j = 0; j < 256; ++j) {
                    const BoundVec3 ray_origin(ray_origin_x, ray_origin_y, ray_origin_z);
                    const FreeVec3 ray_direction(0.0, 0.0, 1.0);
                    const Ray ray(ray_origin, ray_direction);
                    const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
                    assert(!actual_voxels.empty());
                    ray_origin_y = (j != 255) ? ray_origin_y + ray_movement : -1000.0;
                }
                ray_origin_x += ray_movement;
            }
        }
    }

    BENCHMARK(TraversalOne)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalTwo)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalParallelX)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalParallelY)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalParallelZ)->Unit(benchmark::kMillisecond);
    BENCHMARK(MultipleRayNoIntersection)->Unit(benchmark::kMillisecond);
    BENCHMARK(OrthographicRayTracing)->Unit(benchmark::kMillisecond);
    BENCHMARK(OrthographicAllRaysIntersect)->Unit(benchmark::kMillisecond);

} // namespace

BENCHMARK_MAIN();
