% Tests for the PolarCoordinateTraversal function.
% 
% To run: 
%       type <runtests> in the command window.

function tests = polarCoordinateTraversalTest
    tests = functiontests(localfunctions);
end

% Ray doesn't enter circle
function testRayDoesNotEnterCircle(testCase)
    fprintf("Ray doesn't enter circle\n")
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
    fprintf("Ray begins outside and travels through the origin\n")
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

% Ray begins outside of circle and does not go through origin.
function testRayBeginsOutsideCircleAndDoesNotGoThroughOrigin(testCase)
    fprintf("Ray begins outside of circle and does not go through origin\n")
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
        expected_rVoxels     = [1,2,2,1];
        expected_thetaVoxels = [2,2,1,1];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray begins inside of circle and does not go through origin.
function testRayBeginInsideCircleAndDoesNotGoThroughOrigin(testCase)
    fprintf("Ray begins inside of circle and does not go through origin\n")
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
        expected_rVoxels     = [2,3,3,2,2,1];
        expected_thetaVoxels = [2,2,1,1,0,0];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray with zero x-/y-direction.
function testZeroX(testCase)
    fprintf("Ray with zero x-direction\n")
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
    fprintf("Ray with zero y-direction\n")
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
    fprintf("Ray with negative x-/y-direction\n")
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
        expected_rVoxels     = [1,2,2,1,1];
        expected_thetaVoxels = [0,0,3,3,2];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray begins at circle center.
function testBeginAtCircleCenterNegativeXY(testCase)
    fprintf("Ray begins at circle center\n")
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
    fprintf("BeginAtCircleCenterPositiveXY\n")
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
    fprintf("Ray begins at a theta boundary (and goes in orthogonal directions\n")
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
        expected_thetaVoxels = [1,1];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray begins at a radial boundary (and goes in orthogonal directions).
function testBeginAtRadialBoundaryInOrthogonal(testCase)
    fprintf("Ray begins at a radial boundary (and goes in orthogonal directions)\n")
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

% Ray goes along a theta boundary.
function testRayGoesAlongThetaBoundary(testCase)
    fprintf("Ray goes along a theta boundary\n")
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

% Ray intersects a “theta boundary” outside the circle, and ensure it does 
% NOT record this as a theta hit.
function testRayIntersectThetaBoundaryOutsideTheCircle(testCase)
    fprintf("Ray intersects a “theta boundary” outside the circle, and ensure it does NOT record this as a theta hit\n")
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [5.0, 5.0];
    ray_direction = [1.0, 0.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [];
        expected_thetaVoxels = [];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end


% Ray that traverses a single voxel.
function testRayTraverseSingleVoxel(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [12.0, 5.0];
    ray_direction = [-1.0, 1.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1];
        expected_thetaVoxels = [2];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray that ends inside the circle.
function testRayEndsInsideCircle(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [12.0, 5.0];
    ray_direction = [2.0, 1.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 5.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,1];
        expected_thetaVoxels = [2,3];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray with positive x- and negative y- directions.
function testRayWithPositiveXNegativeY(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [10.0, 25.0];
    ray_direction = [2.0, -2.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,2,2,1];
        expected_thetaVoxels = [1,1,0,3,3];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray with negative x- and positive y- directions.
function testRayWithNegativeXPositiveY(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [10.0, 25.0];
    ray_direction = [2.0, -2.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,2,2,1];
        expected_thetaVoxels = [1,1,0,3,3];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray that has a t_begin time different than 0.0
function testRayBeginTimeDiffThan0(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [10.0, 25.0];
    ray_direction = [2.0, -2.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 10.0;
    t_end = 20.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [];
        expected_thetaVoxels = [];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray traversing through circle that has an odd number of angular sections.
function testRayWithOddNumAngularSections(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [5.0, 10.0];
    ray_direction = [3.0, 2.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 3;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,3,3,2,1];
        expected_thetaVoxels = [1,1,1,0,0,0];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Ray that tightly intersects a voxel to determine if floating point tolerances 
% will still deem it an intersection.
function testRaySlightIntersect(testCase)
    min_bound = [0.0, 0.0];
    max_bound = [30.0, 30.0];
    ray_origin = [10.0, 24.96];
    ray_direction = [1.0, 0.0];
    circle_center = [15.0, 15.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,1];
        expected_thetaVoxels = [1,0];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end

% Grid centered at origin; ray traverses through origin
function testRaySlightIntersect(testCase)
    min_bound = [-15, 15];
    max_bound = [15, 15];
    ray_origin = [-13, -13];
    ray_direction = [1.0, 1.0];
    circle_center = [0.0, 0.0];
    circle_max_radius = 10.0;
    num_radial_sections = 3;
    num_angular_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels] = polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose);
        expected_rVoxels     = [1,2,3,3,2,1];
        expected_thetaVoxels = [2,2,2,0,0,0];
        
        verifyEqual(testCase, rVoxels, expected_rVoxels);
        verifyEqual(testCase, thetaVoxels, expected_thetaVoxels);
end
