'''
Testing for cythonized version of the spherical coordinate voxel traversal algorithm.
'''

import unittest
import numpy as np
import cython_SVR

class TestCythonizedSphericalVoxelTraversal(unittest.TestCase):
    # Verifies correctness of the voxel traversal coordinates.
    def verify_voxels(self, voxels, expected_radial_voxels,  expected_theta_voxels,  expected_phi_voxels):
        actual_radial_voxels = []
        actual_theta_voxels = []
        actual_phi_voxels = []
        for x in range(voxels.shape[0]):
            actual_radial_voxels.append(voxels[x, 0])
            actual_theta_voxels.append(voxels[x, 1])
            actual_phi_voxels.append(voxels[x, 2])

        self.assertListEqual(actual_radial_voxels, expected_radial_voxels)
        self.assertListEqual(actual_theta_voxels, expected_theta_voxels)
        self.assertListEqual(actual_phi_voxels, expected_phi_voxels)

    def test_ray_does_not_enter_sphere(self):
        ray_origin = np.array([3.0, 3.0, 3.0])
        ray_direction = np.array([-2.0, -1.3, 1.0])
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([30.0, 30.0, 30.0])
        sphere_center = np.array([15.0, 15.0, 15.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 8
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 15.0
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        assert voxels.size == 0

    def test_sphere_center_at_origin(self):
        ray_origin = np.array([-13.0, -13.0, -13.0])
        ray_direction = np.array([1.0, 1.0, 1.0])
        min_bound = np.array([-20.0, -20.0, -20.0])
        max_bound = np.array([20.0, 20.0, 20.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 4
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0

        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1,2,3,4,4,3,2,1]
        expected_theta_voxels = [2,2,2,2,0,0,0,0]
        expected_phi_voxels = [2,2,2,2,0,0,0,0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_slight_offset_in_XY_plane(self):
        ray_origin = np.array([-13.0, -13.0, -13.0])
        ray_direction = np.array([1.0, 1.5, 1.0])
        min_bound = np.array([-20.0, -20.0, -20.0])
        max_bound = np.array([20.0, 20.0, 20.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 4
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1, 2, 2, 3, 2, 2, 1]
        expected_theta_voxels = [2, 2, 1, 1, 1, 0, 0]
        expected_phi_voxels = [2, 2, 2, 2, 2, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)


    def test_ray_direction_travels_along_X_axis(self):
        ray_origin = np.array([-15.0, 0.0, 0.0])
        ray_direction = np.array([1.0, 0.0, 0.0])
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([30.0, 30.0, 30.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 8
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0

        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1,2,3,4,4,3,2,1]
        expected_theta_voxels = [3,3,3,3,0,0,0,0]
        expected_phi_voxels = [1,1,1,1,0,0,0,0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_direction_travels_along_Y_axis(self):
        ray_origin = np.array([0.0, -15.0, 0.0])
        ray_direction = np.array([0.0, 1.0, 0.0])
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([30.0, 30.0, 30.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 8
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0

        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1,2,3,4,4,3,2,1]
        expected_theta_voxels = [5,5,5,5,1,1,1,1]
        expected_phi_voxels = [0,0,0,0,0,0,0,0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_direction_travels_along_Z_axis(self):
        ray_origin = np.array([0.0, 0.0, -15.0])
        ray_direction = np.array([0.0, 0.0, 1.0])
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([30.0, 30.0, 30.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 8
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0

        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1,2,3,4,4,3,2,1]
        expected_theta_voxels = [0,0,0,0,0,0,0,0]
        expected_phi_voxels = [2,2,2,2,0,0,0,0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_parallel_to_XY_plane(self):
        ray_origin = np.array([-15.0, -15.0, 0.0])
        ray_direction = np.array([1.0, 1.0, 0.0])
        min_bound = np.array([-20.0, -20.0, -20.0])
        max_bound = np.array([20.0, 20.0, 20.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 4
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0

        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1,2,3,4,4,3,2,1]
        expected_theta_voxels = [2,2,2,2,0,0,0,0]
        expected_phi_voxels = [1,1,1,1,0,0,0,0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_parallel_to_XZ_plane(self):
        ray_origin = np.array([-15.0, 0.0, -15.0])
        ray_direction = np.array([1.0, 0.0, 1.0])
        min_bound = np.array([-20.0, -20.0, -20.0])
        max_bound = np.array([20.0, 20.0, 20.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 4
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0

        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1,2,3,4,4,3,2,1]
        expected_theta_voxels = [1,1,1,1,0,0,0,0]
        expected_phi_voxels = [2,2,2,2,0,0,0,0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_parallel_to_YZ_plane(self):
        ray_origin = np.array([0.0, -15.0, -15.0])
        ray_direction = np.array([0.0, 1.0, 1.0])
        min_bound = np.array([-20.0, -20.0, -20.0])
        max_bound = np.array([20.0, 20.0, 20.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 4
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0

        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1,2,3,4,4,3,2,1]
        expected_theta_voxels = [2,2,2,2,0,0,0,0]
        expected_phi_voxels = [2,2,2,2,0,0,0,0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_dir_neg_Y_positive_XZ(self):
        ray_origin = np.array([-13.0, 17.0, -15.0])
        ray_direction = np.array([1.0, -1.2, 1.3])
        min_bound = np.array([-20.0, -20.0, -20.0])
        max_bound = np.array([20.0, 20.0, 20.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 4
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1, 2, 3, 3, 4, 4, 3, 3, 2, 1]
        expected_theta_voxels = [1, 1, 1, 1, 1, 0, 0, 3, 3, 3]
        expected_phi_voxels = [2, 2, 2, 1, 1, 0, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_dir_neg_Z_positive_XY(self):
        ray_origin = np.array([-13.0, -12.0, 15.3])
        ray_direction = np.array([1.4, 2.0, -1.3])
        min_bound = np.array([-20.0, -20.0, -20.0])
        max_bound = np.array([20.0, 20.0, 20.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 4
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1, 1, 2, 2, 1]
        expected_theta_voxels = [2, 1, 1, 0, 0]
        expected_phi_voxels = [1, 1, 1, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_dir_neg_X_positive_YZ(self):
        ray_origin = np.array([13.0, -15.0, -15.0])
        ray_direction = np.array([-1.0, 1.0, 1.0])
        min_bound = np.array([-20.0, -20.0, -20.0])
        max_bound = np.array([20.0, 20.0, 20.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_angular_sections = 4
        num_azimuthal_sections = 4
        t_begin = 0.0
        t_end = 30.0
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                                  num_angular_sections, num_azimuthal_sections, sphere_center,
                                                  sphere_max_radius, t_begin, t_end)
        expected_radial_voxels = [1, 2, 3, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [3, 3, 3, 2, 2, 1, 1, 1, 1]
        expected_phi_voxels = [3, 3, 3, 2, 2, 1, 1, 1, 1]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

if __name__ == '__main__':
    unittest.main()