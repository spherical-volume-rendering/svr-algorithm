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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
        sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, []);
    verifyEqual(testCase, thetaVoxels, []);
    verifyEqual(testCase, phiVoxels, []);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,2,2,0,0,0,0]);
    verifyEqual(testCase, phiVoxels, [2,2,2,2,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, [1,2,2,3,2,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,1,1,1,0,0]);
    verifyEqual(testCase, phiVoxels, [2,2,2,2,2,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
     sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [0,0,0,0,0,0,0,0]);
    %   Note that phi_voxels transition is somewhat ambiguous, as long as
    %   voxels both lie along the z-axis this should be counted as correct
    verifyEqual(testCase, phiVoxels, [2,2,2,2,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
     sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
        verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [3,3,3,3,0,0,0,0]);
    verifyEqual(testCase, phiVoxels,  [1,1,1,1,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
     sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [5,5,5,5,1,1,1,1]);
    verifyEqual(testCase, phiVoxels,  [0,0,0,0,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,2,2,0,0,0,0]);
    verifyEqual(testCase, phiVoxels, [1,1,1,1,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [1,1,1,1,0,0,0,0]);
    verifyEqual(testCase, phiVoxels, [2,2,2,2,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels, [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,2,2,0,0,0,0]);
    verifyEqual(testCase, phiVoxels, [2,2,2,2,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels,     [1,2,3,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [3,3,3,2,2,1,1,1,1]);
    verifyEqual(testCase, phiVoxels,   [3,3,3,2,2,1,1,1,1]);
    tRelError = (tTest-tTraversal)/(tTest^2)
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels,     [1,2,3,3,4,4,3,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [1,1,1,1,1,0,0,3,3,3]);
    verifyEqual(testCase, phiVoxels,   [2,2,2,1,1,0,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
end

% Ray with negative Z positive XY direction
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
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels,     [1,1,2,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,1,1,0,0]);
    verifyEqual(testCase, phiVoxels,   [1,1,1,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
end

% Ray with negative XYZ direction
function testRayNegXYZ(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [15.0, 12.0, 15.0];
    ray_direction = [-1.4, -2.0, -1.3];
    sphere_center = [0.0, 0.0, 0.0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 4;
    num_angular_sections = 4;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels,     [1,1,2,1,1]);
    verifyEqual(testCase, thetaVoxels, [0,3,3,3,2]);
    verifyEqual(testCase, phiVoxels,   [0,0,0,0,1]);
    tRelError = (tTest-tTraversal)/(tTest^2)
end

% Ray with odd number of radial sections
function testOddRadialSections(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-15.0, -15.0, -15.0];
    ray_direction = [1.2, 1.0, 1.3];
    sphere_center = [0.0, 0.0, 0.0];
    sphere_max_radius = 9.0;
    
    num_radial_sections = 3;
    num_angular_sections = 4;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels,     [1,2,2,3,3,2,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,2,2,3,3,0,0]);
    verifyEqual(testCase, phiVoxels,   [2,2,1,1,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
end

% Ray with odd number of angular sections
function testOddAngularSections(testCase)
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
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);

    verifyEqual(testCase, rVoxels,     [1,2,2,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [1,1,1,1,0,0]);
    verifyEqual(testCase, phiVoxels,   [2,2,1,1,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
end

% Ray with odd number of azimuthal sections
function testOddAziSections(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-15.0, -15.0, -15.0];
    ray_direction = [1.0, 1.0, 1.0];
    sphere_center = [0.0, 0.0, 0.0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 4;
    num_angular_sections = 4;
    num_azimuthal_sections = 3;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels,     [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,2,2,0,0,0,0]);
    verifyEqual(testCase, phiVoxels,   [1,1,1,1,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
end

% Ray with large number of radial sections
function testLargeRadialSections(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-15.0, -15.0, -15.0];
    ray_direction = [1.0, 1.0, 1.0];
    sphere_center = [0.0, 0.0, 0.0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 40;
    num_angular_sections = 4;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels,     [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,...
    39,40,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]);
    tRelError = (tTest-tTraversal)/(tTest^2)
end

% Ray with large number of angular sections
function testLargeAngularSections(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-15.0, -15.0, -15.0];
    ray_direction = [1.0, 1.0, 1.0];
    sphere_center = [0.0, 0.0, 0.0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 4;
    num_angular_sections = 40;
    num_azimuthal_sections = 4;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels,     [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [24,24,24,24,4,4,4,4]);
    verifyEqual(testCase, phiVoxels,   [2,2,2,2,0,0,0,0]);

    tRelError = (tTest-tTraversal)/(tTest^2)
end

% Ray with large number of azimuthal sections
function testLargeAziSections(testCase)
    min_bound = [-20, -20.0, -20.0];
    max_bound = [20.0, 20.0, 20.0];
    ray_origin = [-15.0, -15.0, -15.0];
    ray_direction = [1.0, 1.0, 1.0];
    sphere_center = [0.0, 0.0, 0.0];
    sphere_max_radius = 10.0;
    
    num_radial_sections = 4;
    num_angular_sections = 4;
    num_azimuthal_sections = 40;
    t_begin = 0.0;
    t_end = 30.0;
    verbose = false;
    
    [rVoxels, thetaVoxels, phiVoxels, tTest, tTraversal] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, ...
    sphere_center, sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose);
    verifyEqual(testCase, rVoxels,     [1,2,3,4,4,3,2,1]);
    verifyEqual(testCase, thetaVoxels, [2,2,2,2,0,0,0,0]);
    verifyEqual(testCase, phiVoxels,   [24,24,24,24,4,4,4,4]);

    tRelError = (tTest-tTraversal)/(tTest^2)
end