# disutils: language = c++

import numpy as np
cimport numpy as np
cimport cython
from libcpp.vector cimport vector

cdef extern from "../spherical_volume_rendering_util.h" namespace "svr":
    cdef cppclass SphericalVoxel:
        int radial, polar, azimuthal
    cdef cppclass SphereBound:
        double radial, polar, azimuthal

    vector[SphericalVoxel] walkSphericalVolume(double *ray_origin, double *ray_direction,
                                               double *min_bound, double *max_bound,
                                               size_t num_radial_voxels, size_t num_polar_voxels,
                                               size_t num_azimuthal_voxels, double *sphere_center,
                                               double max_t)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def walk_spherical_volume(np.ndarray[np.float64_t, ndim=1, mode="c"] ray_origin,
                          np.ndarray[np.float64_t, ndim=1, mode="c"] ray_direction,
                          np.ndarray[np.float64_t, ndim=1, mode="c"] min_bound,
                          np.ndarray[np.float64_t, ndim=1, mode="c"] max_bound,
                          int num_radial_voxels, int num_polar_voxels, int num_azimuthal_voxels,
                          np.ndarray[np.float64_t, ndim=1, mode="c"] sphere_center, np.float64_t max_t = 1.0):
    '''
    Spherical Coordinate Voxel Traversal Algorithm
    Cythonized version of the Spherical Coordinate Voxel Traversal Algorithm.
    Arguments:
           ray_origin: The 3-dimensional (x,y,z) origin of the ray.
           ray_direction: The 3-dimensional (x,y,z) unit direction of the ray.
           min_bound: The minimum boundary of the sectored sphere in the form (radial, theta, phi).
           max_bound: The maximum boundary of the sectored sphere in the form (radial, theta, phi).
           num_radial_voxels: The number of radial voxels.
           num_polar_voxels: The number of polar voxels.
           num_azimuthal_voxels: The number of azimuthal voxels.
           sphere_center: The 3-dimensional (x,y,z) center of the sphere.
           max_t: The unitized maximum time of ray traversal. Defaulted to 1.0
    Returns:
           A numpy array of the spherical voxel coordinates.
           The voxel coordinates are as follows:
             For coordinate i in numpy array v:
             v[i,0] = radial_voxel
             v[i,1] = polar_voxel
             v[i,2] = azimuthal_voxel
    Notes:
        - If one wants to traverse the entire sphere, the min_bound = [0, 0, 0]
          and max_bound = [SPHERE_MAX_RADIUS, 2*pi, 2*pi]. Similarly, if one wants to traverse
          the upper hemisphere, max_bound = SPHERE_MAX_RADIUS, 2*pi, pi].
        - Code must be compiled before use:
           > python3 cython_SVR_setup.py build_ext --inplace
    '''
    assert(ray_origin.size == 3)
    assert(ray_direction.size == 3)
    assert(sphere_center.size == 3)
    assert(min_bound.size == 3)
    assert(max_bound.size == 3)

    cdef vector[SphericalVoxel] voxels = walkSphericalVolume(&ray_origin[0], &ray_direction[0],
                                                             &min_bound[0], &max_bound[0],
                                                             num_radial_voxels, num_polar_voxels,
                                                             num_azimuthal_voxels, &sphere_center[0],
                                                             max_t)
    cdef np.ndarray cyVoxels = np.empty((voxels.size(), 3), dtype=int)
    for i in range(voxels.size()):
        cyVoxels[i,0] = voxels[i].radial
        cyVoxels[i,1] = voxels[i].polar
        cyVoxels[i,2] = voxels[i].azimuthal
    return cyVoxels