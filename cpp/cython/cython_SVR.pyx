# disutils: language = c++

import numpy as np
cimport numpy as np
cimport cython
from libcpp.vector cimport vector

cdef extern from "../spherical_volume_rendering_util.h" namespace "svr":
    cdef cppclass SphericalVoxel:
        int radial_voxel, angular_voxel, azimuthal_voxel
    cdef cppclass SphereBound:
        double radial, angular, azimuthal

    vector[SphericalVoxel] walkSphericalVolume(double *ray_origin, double *ray_direction,
                                               size_t num_radial_voxels, size_t num_angular_voxels,
                                               size_t num_azimuthal_voxels, double *sphere_center,
                                               double sphere_max_radius, double *min_bound, double *max_bound,
                                               double t_begin, double t_end)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def walk_spherical_volume(np.ndarray[np.float64_t, ndim=1, mode="c"] ray_origin,
                          np.ndarray[np.float64_t, ndim=1, mode="c"] ray_direction,
                          int num_radial_voxels, int num_angular_voxels, int num_azimuthal_voxels,
                          np.ndarray[np.float64_t, ndim=1, mode="c"] sphere_center,
                          np.float64_t sphere_max_radius,
                          np.ndarray[np.float64_t, ndim=1, mode="c"] min_bound,
                          np.ndarray[np.float64_t, ndim=1, mode="c"] max_bound,
                          np.float64_t t_begin, np.float64_t t_end):
    '''
    Spherical Coordinate Voxel Traversal Algorithm
    Cythonized version of the Spherical Coordinate Voxel Traversal Algorithm.
    Arguments:
           ray_origin: The 3-dimensional (x,y,z) origin of the ray.
           ray_direction: The 3-dimensional (x,y,z) direction of the ray.
           num_radial_voxels: The number of radial voxels.
           num_angular_voxels: The number of angular voxels.
           num_azimuthal_voxels: The number of azimuthal voxels.
           sphere_center: The 3-dimensional (x,y,z) center of the sphere.
           min_bound: The minimum boundary of the sectored sphere in the form (radial, theta, phi).
           max_bound: The maximum boundary of the sectored sphere in the form (radial, theta, phi).
           sphere_max_radius: The maximum radius of the sphere.
           t_begin: The beginning time of the ray.
           t_end: The end time of the ray.
    Returns:
           A numpy array of the spherical voxel coordinates.
           The voxel coordinates are as follows:
             For coordinate i in numpy array v:
             v[i,0] = radial_voxel
             v[i,1] = angular_voxel
             v[i,2] = azimuthal_voxel
    Notes:
        - If one wants to traverse the entire sphere, the min_bound = [0, 0, 0]
          and max_bound = [SPHERE_MAX_RADIUS, 2*pi, 2*pi]. Similarly, if one wants to traverse the upper hemisphere,
          max_bound = SPHERE_MAX_RADIUS, 2*pi, pi].
        - Code must be compiled before use:
           > python cython_SVR_setup.py build_ext --inplace
    '''
    assert(ray_origin.size == 3)
    assert(ray_direction.size == 3)
    assert(sphere_center.size == 3)
    assert(min_bound.size == 3)
    assert(max_bound.size == 3)

    cdef vector[SphericalVoxel] voxels = walkSphericalVolume(&ray_origin[0], &ray_direction[0],
                                                             num_radial_voxels, num_angular_voxels,
                                                             num_azimuthal_voxels, &sphere_center[0],
                                                             sphere_max_radius, &min_bound[0], &max_bound[0],
                                                             t_begin, t_end)
    cdef np.ndarray cyVoxels = np.empty((voxels.size(), 3), dtype=int)
    for i in range(voxels.size()):
        cyVoxels[i,0] = voxels[i].radial_voxel
        cyVoxels[i,1] = voxels[i].angular_voxel
        cyVoxels[i,2] = voxels[i].azimuthal_voxel
    return cyVoxels