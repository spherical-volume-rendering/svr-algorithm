% Tests for the PolarCoordinateTraversal function.
% It currently uses polarCoordinateTraversal_old.
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
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal_old(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, []);
    verifyEqual(testCase, thetaVoxels, []);
end

% Ray begins outside and travels through the origin.
function testRayEntersCircleAndGoesThroughOrigin(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [3.0, 5.0];
    ray_direction = [1.0, 1.0];  
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 2;
    num_angular_sections = 4; 
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal_old(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
    expected_rVoxels     = [1,2,2,2,1];
    expected_thetaVoxels = [2,2,1,0,0];
    
    verifyEqual(testCase, rVoxels, expected_rVoxels);
    verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% TODOs:
% Ray begins outside of circle and does not go through origin.
% Ray with zero x-/y-direction.
% Ray with negative x-/y-direction.
% Ray begins at circle center.
% Ray begins at a theta boundary (and goes in orthogonal directions).
% Ray begins at a radial boundary (and does in orthogonal directions).
