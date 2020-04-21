# Fast Voxel Traversal Algorithm Over Spherical Grids
![Example ray tracing in spherical coordinates](images/polar_view_next_to_spherical_image.png)

## About
This project extends the [yt](https://yt-project.org/) open-source data analysis and visualization package, providing an enhanced, integrated user interface for data exploration and enabling the visualization of physical data that comes from non-cartesian grids. Currently, yt implements a fast voxel traversal over a cartesian coordinate grid. The objective is to develop a fast voxel traversal over a spherical coordinate grid, based on ideas from Amanatides and Woo’s seminal paper on fast voxel traversal for ray tracing.

## Authors
- Chris Gyurgyik (cpg49 at cornell.edu)
- Ariel Kellison (ak2485 at cornell.edu)
- Youhan Yuan (yy435 at cornell.edu)

### Algorithm Team Links
- [Fast Voxel Traversal Algorithm Overview](https://docs.google.com/document/d/1QvWw81A0T5vcMAt1WElDeSdBmsw0KJvJdYNr7XfRHfw/edit)
- [Modern C++ implementation of "A Fast Voxel Traversal Algorithm"](https://github.com/cgyurgyik/fast-voxel-traversal-algorithm)

### Project Links
- [Initial Proposal](https://hackmd.io/VRyhXnAFQyaCytWCdKe_1Q)
- [Feasibility Study](https://docs.google.com/document/d/1MbGmy5cSSesI0oUCWHxpiwcHEw6kqd79AV1XZW-rEZo/edit)
- [Progress Report 1](https://docs.google.com/document/d/1ixD7XNu39kwwXhvQooMNb79x18-GsyMPLodzvwC3X-E/edit?ts=5e5d6f45#)

### References
- John Amanatides and Andrew Woo. A fast voxel traversal algorithm for ray tracing. In Eurographics ’87, pages 3–10, 1987.
- James Foley, Andries van Dam, Steven Feiner & John Hughes, "Clipping Lines" in Computer Graphics (3rd Edition) (2013)
- Paul S. Heckbert, editor. Graphics Gems IV.  Academic Press Professional, Inc., USA, 1994.
- Donald. E. Knuth, 1998, Addison-Wesley Longman, Inc., ISBN 0-201-89684-2, Addison-Wesley Professional; 3rd edition.
- Joseph O'Rourke, "Search and  Intersection" in Computational Geometry in C (2nd Edition) (1998)
