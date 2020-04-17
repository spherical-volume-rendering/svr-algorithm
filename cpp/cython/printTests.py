import numpy as np
import sys

# To verify that the C++ and Cythonized versions both work, this script allows for me to type the inputs once and
# print both tests.

cpp_title = "RayDirectionNegativeXPositiveYZ"
python_title = "test_ray_dir_neg_X_positive_YZ"
min_bound = np.array([-20.0, -20.0, -20.0])
max_bound = np.array([20.0, 20.0, 20.0])
sphere_center = np.array([0.0,0.0,0.0])
sphere_max_radius = 10.0
num_radial_sections = 4
num_angular_sections = 4
num_azimuthal_sections = 4
ray_origin = np.array([13.0, -15.0, -15.0])
ray_dir = np.array([-1.0, 1.0, 1.0])
t_begin = 0.0
t_end = 30.0
expected_radial_voxels = str([1,2,3,3,4,4,3,2,1])[1:-1]
expected_angular_voxels = str([3,3,3,2,2,1,1,1,1])[1:-1]
expected_azimuthal_voxels = str([3,3,3,2,2,1,1,1,1])[1:-1]

print("TEST(SphericalCoordinateTraversal, {0}) {{".format(cpp_title))
print("    const BoundVec3 min_bound({0}, {1}, {2});".format(min_bound[0], min_bound[1], min_bound[2]))
print("    const BoundVec3 max_bound({0}, {1}, {2});".format(max_bound[0], max_bound[1], max_bound[2]))
print("    const BoundVec3 sphere_center({0}, {1}, {2});".format(sphere_center[0], sphere_center[1], sphere_center[2]))
print("    const double sphere_max_radius = {0};".format(sphere_max_radius))
print("    const std::size_t num_radial_sections = {0};".format(num_radial_sections))
print("    const std::size_t num_angular_sections = {0};".format(num_angular_sections))
print("    const std::size_t num_azimuthal_sections = {0};".format(num_azimuthal_sections))
print("    const SphericalVoxelGrid grid(min_bound, max_bound, num_radial_sections,")
print("                                  num_angular_sections, num_azimuthal_sections,")
print("                                  sphere_center, sphere_max_radius);")
print("    const BoundVec3 ray_origin({0}, {1}, {2});".format(ray_origin[0], ray_origin[1], ray_origin[2]))
print("    const FreeVec3 ray_direction({0}, {1}, {2});".format(ray_dir[0], ray_dir[1], ray_dir[2]))
print("    const Ray ray(ray_origin, ray_direction);")
print("    const double t_begin = {0};".format(t_begin))
print("    const double t_end = {0};".format(t_end))
print("\n    const auto actual_voxels = sphericalCoordinateVoxelTraversal(ray, grid, t_begin, t_end);")
print("    const std::vector<int> expected_radial_voxels = {{{0}}};".format(expected_radial_voxels))
print("    const std::vector<int> expected_theta_voxels = {{{0}}};".format(expected_angular_voxels))
print("    const std::vector<int> expected_phi_voxels = {{{0}}};".format(expected_azimuthal_voxels))
print("    expectEqualVoxels(actual_voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels);")
print("}")

print("\n")
print("def {0}(self):".format(python_title))
print("    ray_origin = np.array([{0}, {1}, {2}])".format(ray_origin[0], ray_origin[1], ray_origin[2]))
print("    ray_direction = np.array([{0}, {1}, {2}])".format(ray_dir[0], ray_dir[1], ray_dir[2]))
print("    min_bound = np.array([{0}, {1}, {2}])".format(min_bound[0], min_bound[1], min_bound[2]))
print("    max_bound = np.array([{0}, {1}, {2}])".format(max_bound[0], max_bound[1], max_bound[2]))
print("    sphere_center = np.array([{0}, {1}, {2}])".format(sphere_center[0], sphere_center[1], sphere_center[2]))
print("    sphere_max_radius = {0}".format(sphere_max_radius))
print("    num_radial_sections = {0}".format(num_radial_sections))
print("    num_angular_sections = {0}".format(num_angular_sections))
print("    num_azimuthal_sections = {0}".format(num_azimuthal_sections))
print("    t_begin = {0}".format(t_begin))
print("    t_end = {0}".format(t_end))

print("    voxels = CythonSVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,")
print("                                             num_angular_sections, num_azimuthal_sections, sphere_center,")
print("                                             sphere_max_radius, t_begin, t_end)")
print("    expected_radial_voxels = [{0}]".format(expected_radial_voxels))
print("    expected_theta_voxels = [{0}]".format(expected_angular_voxels))
print("    expected_phi_voxels = [{0}]".format(expected_azimuthal_voxels))
print("    self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)")
