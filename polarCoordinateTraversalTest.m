% Tests for the PolarCoordinateTraversal function.
% 
% To run: 
%       type <runtests> in the command window.

function tests = polarCoordinateTraversalTest
    tests = functiontests(localfunctions);
end

% Ray doesn't enter circle
function testRayDoesNotEnterCircle(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [3.0, 3.0];
    ray_direction = [-2.0, -1.3];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 4;
    num_angular_sections = 8;
    t_begin = 0.0;
    t_end = 15.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, []);
    verifyEqual(testCase, thetaVoxels, []);
end

% Ray begins outside and travels through the origin.
function testRayEntersCircleAndGoesThroughOrigin(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [5.0, 5.0];
    ray_direction = [1.0, 1.0];  
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 2;
    num_angular_sections = 4; 
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
    expected_rVoxels     = [1,2,2,1];
    expected_thetaVoxels = [2,2,0,0];
    
    verifyEqual(testCase, rVoxels, expected_rVoxels);
    verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% TODOs:
% Ray begins outside of circle and does not go through origin.
function testRayBeginsOutsideCircleAndDoesNotGoThroughOrigin(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [5.0, 6.0];
    ray_direction = [2.0, 5.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 4;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,2,2,1];
        expected_thetaVoxels = [4,4,3,2,2];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray begins inside of circle and does not go through origin.
function testRayBeginInsideCircleAndDoesNotGoThroughOrigin(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [12.0, 10.0];
    ray_direction = [2.0, 6.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [2,3,3,3,2,1];
        expected_thetaVoxels = [2,2,1,1,0,0];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray with zero x-/y-direction.
function testZeroX(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [10.0, 5.0];
    ray_direction = [0.0, 6.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,2,1];
        expected_thetaVoxels = [2,2,1,1];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

function testZeroY(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [5.0, 10.0];
    ray_direction = [4.0, 0.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,2,1];
        expected_thetaVoxels = [2,2,3,3];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray with negative x-/y-direction.
function testNegativeXandY(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [25.0, 25.0];
    ray_direction = [-4.0, -7.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,2,3,3];
        expected_thetaVoxels = [0,0,3,3,2];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray begins at circle center.
function testBeginAtCircleCenterNegativeXY(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [15.0, 15.0];
    ray_direction = [-4.0, -7.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [3,2,1];
        expected_thetaVoxels = [2,2,2];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

function testBeginAtCircleCenterPositiveXY(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [15.0, 15.0];
    ray_direction = [2.0, 7.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [3,2,1];
        expected_thetaVoxels = [0,0,0];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray begins at a theta boundary (and goes in orthogonal directions).
function testBeginAtThetaBoundaryInOrthogonal(testCase)
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
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [2,1];
        expected_thetaVoxels = [1,1,1];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray begins at a radial boundary (and goes in orthogonal directions).
function testBeginAtRadialBoundaryInOrthogonal(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [22.07, 22.07];
    ray_direction = [-1.0, -1.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,3,3,2,1];
        expected_thetaVoxels = [0,0,0,2,2,2];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Corner Test. Ray goes along a ray boundary
function testRayGoesAlongRayBoundary(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [25.0, 15.0];
    ray_direction = [-1.0, 0.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,3,3,2,1];
        expected_thetaVoxels = [0,0,0,1,1,1];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end