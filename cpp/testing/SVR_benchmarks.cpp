#include <benchmark/benchmark.h>
#include "../spherical_volume_rendering_util.h"

// Benchmarking for the Spherical Volume Rendering algorithm.
// Uses the Google Benchmark library found at:
// https://github.com/google/benchmark

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
            const SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
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
            const SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
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
            const BoundVec3 min_bound4(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound4(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center4(0.0, 0.0, 0.0);
            const double sphere_max_radius4 = 10000.0;
            const std::size_t num_radial_sections4 = 10000;
            const std::size_t num_angular_sections4 = 10000;
            const std::size_t num_azimuthal_sections4 = 10000;
            const SphericalVoxelGrid grid4(min_bound4, max_bound4, num_radial_sections4,
                                           num_angular_sections4,
                                           num_azimuthal_sections4, sphere_center4, sphere_max_radius4);
            const BoundVec3 ray_origin4(-13000.0, 10.0, 10.0);
            const FreeVec3 ray_direction4(1.0, 0.0, 0.0);
            const Ray ray4(ray_origin4, ray_direction4);
            const double t_begin4 = 0.0;
            const double t_end4 = 100000.0;
            const auto actual_voxels4 = sphericalCoordinateVoxelTraversal(ray4, grid4, t_begin4, t_end4);
        }
    }

    static void TraversalParallelY(benchmark::State &state) {
        for (auto _ : state) {
            const BoundVec3 min_bound2(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound2(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center2(0.0, 0.0, 0.0);
            const double sphere_max_radius2 = 10000.0;
            const std::size_t num_radial_sections2 = 10000;
            const std::size_t num_angular_sections2 = 10000;
            const std::size_t num_azimuthal_sections2 = 10000;
            const SphericalVoxelGrid grid2(min_bound2, max_bound2, num_radial_sections2,
                                           num_angular_sections2,
                                           num_azimuthal_sections2, sphere_center2, sphere_max_radius2);
            const BoundVec3 ray_origin2(10.0, -13000.0, 0.0);
            const FreeVec3 ray_direction2(0.0, 1.0, 0.0);
            const Ray ray2(ray_origin2, ray_direction2);
            const double t_begin2 = 0.0;
            const double t_end2 = 100000.0;
            const auto actual_voxels2 = sphericalCoordinateVoxelTraversal(ray2, grid2, t_begin2, t_end2);
        }
    }

    static void TraversalParallelZ(benchmark::State &state) {
        for (auto _ : state) {
            const BoundVec3 min_bound3(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound3(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center3(0.0, 0.0, 0.0);
            const double sphere_max_radius3 = 10000.0;
            const std::size_t num_radial_sections3 = 10000;
            const std::size_t num_angular_sections3 = 10000;
            const std::size_t num_azimuthal_sections3 = 10000;
            const SphericalVoxelGrid grid3(min_bound3, max_bound3, num_radial_sections3,
                                           num_angular_sections3,
                                           num_azimuthal_sections3, sphere_center3, sphere_max_radius3);
            const BoundVec3 ray_origin3(10.0, 0.0, -13000.0);
            const FreeVec3 ray_direction3(0.0, 0.0, 1.0);
            const Ray ray3(ray_origin3, ray_direction3);
            const double t_begin3 = 0.0;
            const double t_end3 = 100000.0;
            const auto actual_voxels3 = sphericalCoordinateVoxelTraversal(ray3, grid3, t_begin3, t_end3);
        }
    }

    static void MultipleRayNoIntersection(benchmark::State &state) {
        for (auto _ : state) {
            const int number_of_runs = 50000;
            const BoundVec3 min_bound5(0.0, 0.0, 0.0);
            const BoundVec3 max_bound5(30.0, 30.0, 30.0);
            const BoundVec3 sphere_center5(15.0, 15.0, 15.0);
            const double sphere_max_radius5 = 10.0;
            const std::size_t num_radial_sections5 = 4;
            const std::size_t num_angular_sections5 = 8;
            const std::size_t num_azimuthal_sections5 = 4;
            const SphericalVoxelGrid grid5(min_bound5, max_bound5, num_radial_sections5,
                                           num_angular_sections5,
                                           num_azimuthal_sections5, sphere_center5, sphere_max_radius5);
            const BoundVec3 ray_origin5(3.0, 3.0, 3.0);
            const FreeVec3 ray_direction5(-2.0, -1.3, 1.0);
            const Ray ray5(ray_origin5, ray_direction5);
            const double t_begin5 = 0.0;
            const double t_end5 = 15.0;
            for (int i = 0; i < number_of_runs; ++i) {
                sphericalCoordinateVoxelTraversal(ray5, grid5, t_begin5, t_end5);
            }
        }
    }

    static void MultipleRayIntersection(benchmark::State &state) {
        for (auto _ : state) {
            const int number_of_rays = 100;
            const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center(0.0, 0.0, 0.0);
            const double sphere_max_radius = 10000.0;
            const std::size_t num_radial_sections = 10000;
            const std::size_t num_angular_sections = 10000;
            const std::size_t num_azimuthal_sections = 10000;
            const SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                          num_angular_sections,
                                          num_azimuthal_sections, sphere_center, sphere_max_radius);
            const double t_begin = 0.0;
            const double t_end = 100000.0;

            double ray_origin_x = -5000.0;
            for (int i = 0; i < number_of_rays; ++i) {
                const BoundVec3 ray_origin(ray_origin_x, 0.0, 0.0);
                const FreeVec3 ray_direction(1.0, 0.0, 0.0);
                const Ray ray(ray_origin, ray_direction);
                const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
                ray_origin_x += 100.0;
            }
        }
    }

    static void OrthographicRayTrace(benchmark::State &state) {
        for (auto _ : state) {
            const BoundVec3 min_bound(-20000.0, -20000.0, -20000.0);
            const BoundVec3 max_bound(20000.0, 20000.0, 20000.0);
            const BoundVec3 sphere_center(0.0, 0.0, 0.0);
            const double sphere_max_radius = 10000.0;
            const std::size_t num_radial_sections = 10000;
            const std::size_t num_angular_sections = 10000;
            const std::size_t num_azimuthal_sections = 10000;
            const SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,
                                          num_angular_sections,
                                          num_azimuthal_sections, sphere_center, sphere_max_radius);
            const double t_begin = 0.0;
            const double t_end = 100000.0;

            double ray_origin_x = -5000.0;
            double ray_origin_y = -5000.0;
            for (int i = 0; i < 100; ++i) {
                for (int j = 0; j < 100; ++j) {
                    const BoundVec3 ray_origin(ray_origin_x, ray_origin_y, 0.0);
                    const FreeVec3 ray_direction(1.0, 0.0, 0.0);
                    const Ray ray(ray_origin, ray_direction);
                    const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);
                    ray_origin_y += 100.0;
                }
                ray_origin_x += 100.0;
            }
        }
    }

    BENCHMARK(TraversalOne)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalTwo)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalParallelX)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalParallelY)->Unit(benchmark::kMillisecond);
    BENCHMARK(TraversalParallelZ)->Unit(benchmark::kMillisecond);
    BENCHMARK(MultipleRayNoIntersection)->Unit(benchmark::kMillisecond);
    BENCHMARK(MultipleRayIntersection)->Unit(benchmark::kMillisecond);
    BENCHMARK(OrthographicRayTrace)->Unit(benchmark::kMillisecond);
}

BENCHMARK_MAIN();
