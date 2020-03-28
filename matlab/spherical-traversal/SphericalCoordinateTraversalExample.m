    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-13.0, -13.0, -13.0];
    ray_direction = [1.0, 1.0, 1.0];
    sphere_center = [0.0, 0.0, 0.0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 4;
    num_angular_sections = 4;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
