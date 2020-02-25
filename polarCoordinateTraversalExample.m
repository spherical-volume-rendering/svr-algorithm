function polarCoordinateTraversalExample
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    
    ray_origin = [7.0, 10.0];
    ray_direction = [1.0, 1.0];
    
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 4;
    num_angular_sections = 8;
    
    t_begin = 0.0;
    t_end = 15.0;
    
    verbose = true;
   
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
end
