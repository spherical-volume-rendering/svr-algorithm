from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

"""
This script creates three plots: the 3D spherical surface, the XY plane, and
the XZ plane; each with the relevant ray or ray projection.
resources:
1. https://stackoverflow.com/questions/42281966/how-to-plot-vectors-in-python-using-matplotlib
2. https://www.geeksforgeeks.org/vector-projection-using-python/
3. https://pundit.pratt.duke.edu/wiki/Python:Plotting_Surfaces

Requires:
   python3
   matplotlib version == 2.0.0
"""
####################################
# BEGIN EDITABLE
# max radius of the sphere
r = 1
# number of radial sections
num_rad = 10
dr = r/num_rad
# number of angular sections
num_ang = 8
dt = 2*np.pi/num_ang
# number of azimuthal sections
num_azi = 8
dp = 2*np.pi/num_azi
# sphere center
origin_sphere = np.array([0,0,0])
# ray Start
origin_ray = np.array([-1,-1,-1])
# ray direction
ray_dir = np.array([1,0.5,1])
# vector n_xy: n_xy is orthogonal to xy plane
n_xy = np.array([0, 0, 1])
# vector n_xz: n_xz is orthogonal to xz plane
n_xz = np.array([0, 1, 0])
# END EDITABLE
####################################

# Sphere plot axes
fig1 = plt.figure(1)
ax1 = fig1.gca(projection='3d')
ax1.set_aspect("equal")

# draw sphere
u, v = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
x = r * np.cos(u)*np.sin(v) + origin_sphere[0]
y = r * np.sin(u)*np.sin(v) + origin_sphere[1]
z = r * np.cos(v) + origin_sphere[2]
# alpha is transparency
ax1.plot_surface(x, y, z, color="b", alpha = 0.5)

# draw the ray origin point
ax1.scatter([origin_ray[0]], [origin_ray[1]], [origin_ray[2]], color="r", s=100)

# plot the XY plane
min_bound = ax1.get_ylim()[0]
max_bound = ax1.get_ylim()[1]
xx, yy = np.meshgrid(np.arange(min_bound, max_bound, .1), np.arange(min_bound, max_bound, .1))
ax1.plot_surface(xx, yy, origin_sphere[2], alpha=0.5)

# draw the ray
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)
# estimate ray end for grid
end_ray = origin_ray + abs(max_bound) * ray_dir
a = Arrow3D([origin_ray[0], end_ray[0]], [origin_ray[1], end_ray[1]], [origin_ray[2], end_ray[2]], mutation_scale=20,
            lw=1, arrowstyle="-|>", color="k")
ax1.add_artist(a)

# XY plane plot axes
fig2 = plt.figure(2)
ax2 = fig2.gca()
ax2.set_aspect("equal")

# draw sphere projection and angular sections
rad = r
while rad > dr:
    x = rad * np.cos(u) + origin_sphere[0]
    y = rad * np.sin(u) + origin_sphere[1]
    ax2.plot(x, y, color="b")
    rad = rad - dr
theta = np.arange(0,2*np.pi,dt)
xpts = r * np.cos(theta) + origin_sphere[0]
ypts = r * np.sin(theta) + origin_sphere[1]
for i in range(len(xpts)):
    ax2.plot([xpts[i],origin_sphere[0]],[ypts[i],origin_sphere[1]],color="b")

# Euclidean norm of n
n_norm = np.sqrt(sum(n_xy**2))

rayProj_n = (np.dot(ray_dir, n_xy)/n_norm**2)*n_xy

# projection of ray on xy plane
projP = ray_dir - rayProj_n
ax2.quiver([origin_ray[0]], [origin_ray[1]], [projP[0]], [projP[1]], color=['k'], angles='xy', scale_units='xy', scale=1)

# XZ plane plot axes
fig3 = plt.figure(3)
ax3 = fig3.gca()
ax3.set_aspect("equal")

# draw sphere projection and angular sections
rad = r
while rad > dr:
    x = rad * np.cos(u) + origin_sphere[0]
    z = rad * np.sin(u) + origin_sphere[2]
    ax3.plot(x, z, color="g")
    rad = rad - dr
phi = np.arange(0,2*np.pi,dp)
xpts = r * np.cos(phi) + origin_sphere[0]
zpts = r * np.sin(phi) + origin_sphere[2]
for i in range(len(xpts)):
    ax3.plot([xpts[i],origin_sphere[0]],[zpts[i],origin_sphere[2]],color="g")

# Euclidean norm of n
n_norm = np.sqrt(sum(n_xz**2))

rayProj_n = (np.dot(ray_dir, n_xz)/n_norm**2)*n_xz

# projection of ray on xy plane
projP = ray_dir - rayProj_n
ax3.quiver([origin_ray[0]], [origin_ray[2]], [projP[0]], [projP[2]], color=['k'], angles='xy', scale_units='xy', scale=1)

plt.show()
