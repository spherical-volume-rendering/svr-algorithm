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
function testRayOffsetInXYPlane(testCase)
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
    verifyEqual(testCase, rVoxels, [1,2,2,3,2,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,1,1,1,0,0]);
    verifyEqual(testCase, phiVoxels, [2,2,2,2,2,0,0]);
end

% Ray travels up z-axis
function testRayTravelsAlongZAxis(testCase)
    min_bound = [0.0, 0.0, 0.0];
    max_bound = [30.0, 30.0, 30.0];
    ray_origin = [0.0, 0.0, -15.0];
    ray_direction = [0, 0, 1.0];
    sphere_center = [0, 0, 0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 4;
    num_angular_sections = 8;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
     sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [0,0,0,0,0,0,0,0]);
    %   Note that phi_voxels transition is somewhat ambiguous, as long as
    %   voxels both lie along the z-axis this should be counted as correct
    verifyEqual(testCase, phiVoxels, [2,2,2,2,0,0,0,0]);
end

% Ray travels along x-axis
function testRayTravelsAlongXAxis(testCase)
    min_bound = [0.0, 0.0, 0.0];
    max_bound = [30.0, 30.0, 30.0];
    ray_origin = [-15.0, 0.0, 0.0];
    ray_direction = [1.0, 0.0, 0.0];
    sphere_center = [0, 0, 0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 4;
    num_angular_sections = 8;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
     sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
        verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [3,3,3,3,0,0,0,0]);
    verifyEqual(testCase, phiVoxels,  [1,1,1,1,0,0,0,0]);
end

% Ray travels along y-axis
function testRayTravelsAlongYAxis(testCase)
    min_bound = [0.0, 0.0, 0.0];
    max_bound = [30.0, 30.0, 30.0];
    ray_origin = [0.0, -15.0, 0.0];
    ray_direction = [0.0, 1.0, 0.0];
    sphere_center = [0, 0, 0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 4;
    num_angular_sections = 8;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
     sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [5,5,5,5,1,1,1,1]);
    verifyEqual(testCase, phiVoxels,  [0,0,0,0,0,0,0,0]);
end

% Ray parallel to the XY Plane.
function testRayParallelToXYPlane(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-15.0, -15.0, 0.0];
    ray_direction = [1.0, 1.0, 0.0];
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
    verifyEqual(testCase, phiVoxels, [1,1,1,1,0,0,0,0]);
end

% Ray parallel to the XZ Plane.
function testRayParallelToXZPlane(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-15.0, 0.0, -15.0];
    ray_direction = [1.0, 0.0, 1.0];
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
    verifyEqual(testCase, thetaVoxels, [1,1,1,1,0,0,0,0]);
    verifyEqual(testCase, phiVoxels, [2,2,2,2,0,0,0,0]);
end

% Ray parallel to the YZ Plane.
function testRayParallelToYZPlane(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [0.0, -15.0, -15.0];
    ray_direction = [0.0, 1.0, 1.0];
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

% Ray with negative X positive YZ direction
function testRayNegXPosYZ(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [13.0, -15.0, -15.0];
    ray_direction = [-1.0, 1.0, 1.0];
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
    verifyEqual(testCase, rVoxels,     [1,2,3,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [3,3,3,2,2,1,1,1,1]);
    verifyEqual(testCase, phiVoxels,   [3,3,3,2,2,1,1,1,1]);
end

% Ray with negative Y positive XZ direction
function testRayNegYPosXZ(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-13.0, 17.0, -15.0];
    ray_direction = [1.0, -1.2, 1.3];
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
    verifyEqual(testCase, rVoxels,     [1,2,3,3,4,4,3,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [1,1,1,1,1,0,0,3,3,3]);
    verifyEqual(testCase, phiVoxels,   [2,2,2,1,1,0,0,0,0,0]);
end

% Ray with negative Z positive XY direction, NOT COMPLETED
function testRayNegZPosXY(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-13.0, -12.0, 15.3];
    ray_direction = [1.4, 2.0, -1.3];
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
    verifyEqual(testCase, rVoxels,     [1,1,2,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,1,1,0,0]);
    verifyEqual(testCase, phiVoxels,   [1,1,1,0,0]);
end

% Ray with negative XY positive Z direction

% Ray with negative YZ positive X direction

% Ray with negative XZ positive Y direction

% Ray with negative XYZ direction

% Ray with odd number of angular sections

% Ray with odd number of azimuthal sections
