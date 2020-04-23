    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-15.0, -15.0, -15.0];
    ray_direction = [1.0, 1.0, 1.3];
    sphere_center = [0.0, 0.0, 0.0];
    sphere_max_radius = 9.0;
    
    num_radial_sections = 4;
    num_angular_sections = 3;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = true;
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose)
   % verifyEqual(testCase, rVoxels,     [1,2,3,3,2,1]);
   % verifyEqual(testCase, thetaVoxels, [1,1,1,0,0,0]);
   % verifyEqual(testCase, phiVoxels,   [2,2,1,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)