'''
Testing for cythonized version of the spherical coorinate voxel traversal algorithm.
'''

import numpy as np
import CythonSVR


def test_ray_does_not_enter_sphere():
    print("Testing [Ray does not enter_sphere] ...")
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
    tol = 10e-16
    voxels = CythonSVR.walk_spherical_volume(ray_origin, ray_direction, min_bound, max_bound, num_radial_sections,
                                 num_angular_sections, num_azimuthal_sections, sphere_center, sphere_max_radius,
                                 t_begin, t_end, tol)
    empty = np.array([])
    assert np.array_equal(voxels, empty)
    print("... [Ray does not enter sphere] passed.")

test_ray_does_not_enter_sphere()