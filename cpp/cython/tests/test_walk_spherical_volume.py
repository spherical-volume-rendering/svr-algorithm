'''
Testing for cythonized version of the spherical coordinate voxel traversal algorithm.
To compile SVR code:
    python3 cython_SVR_setup.py build_ext --inplace
'''

import unittest
import numpy as np
import cython_SVR

class TestWalkSphericalVolume(unittest.TestCase):
    # Verifies correctness of the voxel traversal coordinates.
    def verify_voxels(self, voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels):
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
        sphere_center = np.array([15.0, 15.0, 15.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 8
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        assert voxels.size == 0

    def test_ray_does_not_enter_sphere_tangential_hit(self):
        ray_origin = np.array([-10.0, -10.0, 0.0])
        ray_direction = np.array([0.0, 1.0, 0.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 8
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        assert voxels.size == 0

    def test_sphere_center_at_origin(self):
        ray_origin = np.array([-13.0, -13.0, -13.0])
        ray_direction = np.array([1.0, 1.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        expected_phi_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_sphere_center_not_at_origin(self):
        ray_origin = np.array([-11.0, -11.0, -11.0])
        ray_direction = np.array([1.0, 1.0, 1.0])
        sphere_center = np.array([2.0, 2.0, 2.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        expected_phi_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_within_sphere(self):
        ray_origin = np.array([-3.0, 4.0, 5.0])
        ray_direction = np.array([1.0, -1.0, -1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [2, 3, 4, 4, 4, 4, 3, 2, 1]
        expected_theta_voxels = [1, 1, 1, 0, 3, 3, 3, 3, 3]
        expected_phi_voxels = [1, 1, 1, 0, 0, 3, 3, 3, 3]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_within_sphere_and_begin_time_is_not_zero(self):
        ray_origin = np.array([-3.0, 4.0, 5.0])
        ray_direction = np.array([1.0, -1.0, -1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4
        t_begin = 5.0
        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [4, 3, 2, 1]
        expected_theta_voxels = [3, 3, 3, 3]
        expected_phi_voxels = [0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_ends_within_sphere(self):
        ray_origin = np.array([13.0, -15.0, 16.0])
        ray_direction = np.array([-1.5, 1.2, -1.5])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 2, 3]
        expected_theta_voxels = [3, 3, 2, 2]
        expected_phi_voxels = [0, 0, 1, 1]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_and_ends_within_sphere(self):
        ray_origin = np.array([-3.0, 4.0, 5.0])
        ray_direction = np.array([1.0, -1.0, -1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [2, 3, 4, 4, 4]
        expected_theta_voxels = [1, 1, 1, 0, 3]
        expected_phi_voxels = [1, 1, 1, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_and_ends_within_sphere_not_centered_at_origin(self):
        ray_origin = np.array([-1.0, 7.0, 7.0])
        ray_direction = np.array([1.0, -1.0, -1.0])
        sphere_center = np.array([2.0, 3.0, 2.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [2, 3, 4, 4, 4]
        expected_theta_voxels = [1, 1, 1, 0, 3]
        expected_phi_voxels = [1, 1, 1, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_slight_offset_in_XY_plane(self):
        ray_origin = np.array([-13.0, -13.0, -13.0])
        ray_direction = np.array([1.0, 1.5, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 2, 3, 2, 2, 1]
        expected_theta_voxels = [2, 2, 1, 1, 1, 0, 0]
        expected_phi_voxels = [2, 2, 2, 2, 2, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_direction_travels_along_X_axis(self):
        ray_origin = np.array([-15.0, 0.0, 0.0])
        ray_direction = np.array([1.0, 0.0, 0.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 8
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [3, 3, 3, 3, 0, 0, 0, 0]
        expected_phi_voxels = [1, 1, 1, 1, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_direction_travels_along_Y_axis(self):
        ray_origin = np.array([0.0, -15.0, 0.0])
        ray_direction = np.array([0.0, 1.0, 0.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 8
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [5, 5, 5, 5, 1, 1, 1, 1]
        expected_phi_voxels = [0, 0, 0, 0, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_direction_travels_along_Z_axis(self):
        ray_origin = np.array([0.0, 0.0, -15.0])
        ray_direction = np.array([0.0, 0.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 8
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [0, 0, 0, 0, 0, 0, 0, 0]
        expected_phi_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_parallel_to_XY_plane(self):
        ray_origin = np.array([-15.0, -15.0, 0.0])
        ray_direction = np.array([1.0, 1.0, 0.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        expected_phi_voxels = [1, 1, 1, 1, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_parallel_to_XZ_plane(self):
        ray_origin = np.array([-15.0, 0.0, -15.0])
        ray_direction = np.array([1.0, 0.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [1, 1, 1, 1, 0, 0, 0, 0]
        expected_phi_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_parallel_to_YZ_plane(self):
        ray_origin = np.array([0.0, -15.0, -15.0])
        ray_direction = np.array([0.0, 1.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        expected_phi_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_dir_neg_Y_positive_XZ(self):
        ray_origin = np.array([-13.0, 17.0, -15.0])
        ray_direction = np.array([1.0, -1.2, 1.3])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 3, 4, 4, 3, 3, 2, 1]
        expected_theta_voxels = [1, 1, 1, 1, 1, 0, 0, 3, 3, 3]
        expected_phi_voxels = [2, 2, 2, 1, 1, 0, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_dir_neg_Z_positive_XY(self):
        ray_origin = np.array([-13.0, -12.0, 15.3])
        ray_direction = np.array([1.4, 2.0, -1.3])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 1, 2, 2, 1]
        expected_theta_voxels = [2, 1, 1, 0, 0]
        expected_phi_voxels = [1, 1, 1, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_dir_neg_X_positive_YZ(self):
        ray_origin = np.array([13.0, -15.0, -15.0])
        ray_direction = np.array([-1.0, 1.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [3, 3, 3, 2, 2, 1, 1, 1, 1]
        expected_phi_voxels = [3, 3, 3, 2, 2, 1, 1, 1, 1]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_dir_neg_XYZ(self):
        ray_origin = np.array([15.0, 12.0, 15.0])
        ray_direction = np.array([-1.4, -2.0, -1.3])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 1, 2, 1, 1]
        expected_theta_voxels = [0, 3, 3, 3, 2]
        expected_phi_voxels = [0, 0, 0, 0, 1]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_odd_number_angular_sections(self):
        ray_origin = np.array([-15.0, -15.0, -15.0])
        ray_direction = np.array([1.0, 1.0, 1.3])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 9.0
        num_radial_sections = 4
        num_polar_sections = 3
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 2, 3, 2, 1]
        expected_theta_voxels = [1, 1, 1, 1, 0, 0]
        expected_phi_voxels = [2, 2, 1, 1, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_odd_number_azimuthal_sections(self):
        ray_origin = np.array([-15.0, -15.0, -15.0])
        ray_direction = np.array([1.0, 1.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 3

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        expected_phi_voxels = [1, 1, 1, 1, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_large_number_of_radial_sections(self):
        ray_origin = np.array([-15.0, -15.0, -15.0])
        ray_direction = np.array([1.0, 1.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 40
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                                  18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
                                  33, 34, 35, 36, 37, 38, 39, 40, 40, 39, 38, 37, 36, 35, 34,
                                  33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19,
                                  18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        expected_theta_voxels = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        expected_phi_voxels = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                               2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_large_number_of_angular_sections(self):
        ray_origin = np.array([-15.0, -15.0, -15.0])
        ray_direction = np.array([1.0, 1.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 40
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [24, 24, 24, 24, 4, 4, 4, 4]
        expected_phi_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_large_number_of_azimuthal_sections(self):
        ray_origin = np.array([-15.0, -15.0, -15.0])
        ray_direction = np.array([1.0, 1.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 40

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [2, 2, 2, 2, 0, 0, 0, 0]
        expected_phi_voxels = [24, 24, 24, 24, 4, 4, 4, 4]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_time_begins_is_not_zero(self):
        ray_origin = np.array([-15.0, 15.0, 15.0])
        ray_direction = np.array([1.0, -1.0, -1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4
        1
        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 4, 4, 3, 2, 1]
        expected_theta_voxels = [1, 1, 1, 1, 3, 3, 3, 3]
        expected_phi_voxels = [1, 1, 1, 1, 3, 3, 3, 3]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_in_outermost_radius_and_ends_within_sphere(self):
        ray_origin = np.array([-4.0, -4.0, -6.0])
        ray_direction = np.array([1.3, 1.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 3, 4, 4]
        expected_theta_voxels = [2, 2, 2, 3, 3, 0]
        expected_phi_voxels = [2, 2, 2, 3, 3, 3]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_at_sphere_origin(self):
        ray_origin = np.array([0.0, 0.0, 0.0])
        ray_direction = np.array([-1.5, 1.2, -1.5])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [4, 3, 2, 1]
        expected_theta_voxels = [1, 1, 1, 1]
        expected_phi_voxels = [2, 2, 2, 2]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_past_sphere_origin_one(self):
        ray_origin = np.array([-3.0, 2.4, -3.0])
        ray_direction = np.array([-1.5, 1.2, -1.5])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [3, 2, 1]
        expected_theta_voxels = [1, 1, 1]
        expected_phi_voxels = [2, 2, 2]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_past_sphere_origin_two(self):
        ray_origin = np.array([-4.5, 3.6, -4.5])
        ray_direction = np.array([-1.5, 1.2, -1.5])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [2, 1]
        expected_theta_voxels = [1, 1]
        expected_phi_voxels = [2, 2]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_past_sphere_origin_three(self):
        ray_origin = np.array([-6.0, 4.8, -6.0])
        ray_direction = np.array([-1.5, 1.2, -1.5])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1]
        expected_theta_voxels = [1]
        expected_phi_voxels = [2]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_begins_past_sphere_origin_four(self):
        ray_origin = np.array([-7.5, 6.0, -7.5])
        ray_direction = np.array([-1.5, 1.2, -1.5])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = []
        expected_theta_voxels = []
        expected_phi_voxels = []
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_tangential_hit(self):
        ray_origin = np.array([-5.0, 0.0, 10.0])
        ray_direction = np.array([0.0, 0.0, -1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 2, 1]
        expected_theta_voxels = [1, 1, 1, 1]
        expected_phi_voxels = [1, 1, 2, 2]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_ray_tangential_hit_two(self):
        ray_origin = np.array([-2.5, 0.0, 10.0])
        ray_direction = np.array([0.0, 0.0, -1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 3, 2, 1]
        expected_theta_voxels = [1, 1, 1, 1, 1, 1]
        expected_phi_voxels = [1, 1, 1, 2, 2, 2]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_tangential_hit_no_double_intersection_with_same_voxel(self):
        ray_origin = np.array([-2.5, 0.0, 10.0])
        ray_direction = np.array([0.0, 0.0, -1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 1
        num_azimuthal_sections = 1

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 2, 1]
        expected_theta_voxels = [0, 0, 0, 0, 0]
        expected_phi_voxels = [0, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_nearly_tangential_hit(self):
        ray_origin = np.array([-5.01, 0.0, 10.0])
        ray_direction = np.array([0.0, 0.0, -1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 4
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, 2 * np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 2, 1]
        expected_theta_voxels = [1, 1, 1, 1]
        expected_phi_voxels = [1, 1, 2, 2]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_upper_hemisphere_hit(self):
        ray_origin = np.array([-11.0, 2.0, 1.0])
        ray_direction = np.array([1.0, 0.0, 0.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 8
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        expected_radial_voxels = [1, 2, 3, 3, 4, 4, 4, 4, 3, 3, 2, 1]
        expected_theta_voxels = [3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0]
        expected_phi_voxels = [3, 3, 3, 3, 3, 2, 1, 0, 0, 0, 0, 0]
        self.verify_voxels(voxels, expected_radial_voxels, expected_theta_voxels, expected_phi_voxels)

    def test_upper_hemisphere_miss(self):
        ray_origin = np.array([-5.0, -5.0, -5.0])
        ray_direction = np.array([1.0, 0.0, 0.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10.0
        num_radial_sections = 4
        num_polar_sections = 8
        num_azimuthal_sections = 4

        t_end = 1.0
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        assert voxels.size == 0

    def test_avoid_ray_stepping_to_radial_voxel_zero(self):
        ray_origin = np.array([-984.375, 250.0, -10001.0])
        ray_direction = np.array([0.0, 0.0, 1.0])
        sphere_center = np.array([0.0, 0.0, 0.0])
        sphere_max_radius = 10e3
        num_radial_sections = 128
        num_polar_sections = 128
        num_azimuthal_sections = 128

        t_end = sphere_max_radius * 3
        min_bound = np.array([0.0, 0.0, 0.0])
        max_bound = np.array([sphere_max_radius, 2 * np.pi, np.pi])
        voxels = cython_SVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound,
                                                  num_radial_sections, num_polar_sections, num_azimuthal_sections,
                                                  sphere_center, t_end)
        last_radial_voxel = voxels[voxels[0].size - 1][0]
        assert (last_radial_voxel != 0)


if __name__ == '__main__':
    unittest.main()
