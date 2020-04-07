% Tests for the SphericalCoordinateTraversal function.
% 
% To run: 
%       type <runtests> in the command window.

function tests = sphericalCoordinateTraversalTest
    tests = functiontests(localfunctions);
end

% Ray doesn't enter sphere
function testRayDoesNotEnterSphere(testCase)
    fprintf("Ray doesn't enter sphere\n")
    min_bound = [0.0, 0.0, 0.0];
    max_bound = [30.0, 30.0, 30.0];
    ray_origin = [3.0, 3.0, 3.0];
    ray_direction = [-2.0, -1.3, 1.0];
    sphere_center = [15.0, 15.0, 15.0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 4;
    num_angular_sections = 8;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 15.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, []);
    verifyEqual(testCase, thetaVoxels, []);
    verifyEqual(testCase, phiVoxels, []);
end

% Ray passes through origin
function testSphereCenteredAtOrigin(testCase)
    fprintf("Sphere at Origin\n")
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
    verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,2,2,0,0,0,0]);
    verifyEqual(testCase, phiVoxels, [2,2,2,2,0,0,0,0]);
end

% Change direction of ray slightly in the XY plane.
function testRayOffsetinXYPlane(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-13.0, -13.0, -13.0];
    ray_direction = [1.0, 1.5, 1.0];
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
    verifyEqual(testCase, rVoxels, [1,2,2,3,2,2,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,1,1,1,1,0,0]);
    verifyEqual(testCase, phiVoxels, [2,2,2,2,2,0,0,0]);
end
