    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [10.0, 15.0];
    ray_direction = [0.0, 7.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = true;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose)
        expected_rVoxels     = [2,1]
        expected_thetaVoxels = [1,1]