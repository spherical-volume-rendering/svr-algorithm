function [radial_voxels, angular_voxels, azimuthal_voxels, time_test, traversal_time] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, sphere_center, ...
    sphere_max_radius, num_radial_sections, num_polar_sections, num_azimuthal_sections, t_begin, t_end, verbose)
% Input:
%    min_bound: The lower left corner of the bounding box.
%    max_bound: The upper right corner of the bounding box.
%    ray origin: The origin of the ray in (x, y, z) coordinates.
%    ray direction: The direction of the ray in (x, y, z) coordinates.
%    sphere_center: The x, y location of the center of the sphere.
%    sphere_max_radius: The largest that encompasses the sphere.
%    num_radial_sections: The number of radial sections in the sphere.
%    num_polar_sections: The number of angular sections in the sphere.
%                          The angular sections lie on the x-y plane.
%    num_azimuthal_sections: The number of azimuthal sections in the sphere.
%                            The azimuthal sections lie on the x-z plane.
%    t_begin: The beginning time of the ray.
%    t_end: The end time of the ray.
%
% Requires:
%    max_bound > min_bound
%    The entire sphere is within min_bound and max_bound.
%    t_end > t_begin >= 0.0
%    sphere_max_radius > 0
%    num_radial_sections > 0
%    num_polar_sections > 0
%    num_azimuthal_sections > 0
%
% Returns:
%    radial_voxels: A list of the radial voxels that were hit by the ray.
%    angular_voxels: A list of the angular voxels that were hit by the ray.
%    azimuthal_voxels: A list of the phi voxels that were hit by the ray.
%    time_test: A float of expected traversal time of the ray in the grid;
%    traversal_test: A float of actual total time spent on traversal in the grid.
%
%    Note: These lists, used in conjunction, will produce the path of the ray
%    through the voxels using each voxel transition. For example,
%    [radial_voxels(1), angular_voxels(1), azimuthal_voxels(1)]
%    is the first voxel the ray travels through. If the next traversal is a radial hit,
%    angular_voxels(2) and azimuthal_voxels(2) will remain the same.

angular_voxels = [];
radial_voxels = [];
azimuthal_voxels = [];

% INITIALIZATION PHASE
% Determine if the ray intersects the grid; Find the intersection times 
% for the ray and the outermost radial shell.
[t1, t2] = radial_intersection_times(ray_origin, ray_direction,...
        sphere_center, sphere_max_radius);

% The ray may not intersect the grid at all.
% In particular, if the ray is outside the grid at t_begin.
if t1 < t_begin && t2 < t_begin 
    time_test = 1.0;
    traversal_time = 1.0;
    if verbose
        fprintf("\nRay does not intersect spherical grid for t_begin.")
    end
    return;
end

% It may be a tangent hit.
 if approximatelyEqual(t1,t2,1e-12,1e-8)
    time_test = 1.0;
    traversal_time = 1.0;
    if verbose
        fprintf("\nTangent hit.")
    end
    return;
end
% Calculate initial radial voxel ID and radius of initial voxel.
[current_voxel_ID_r, r, delta_radius] = radial_init(ray_origin, ray_direction,...
    sphere_center, sphere_max_radius, num_radial_sections, t_begin);

% Calculate initial angular voxel ID and the initial point of intersection
% between the 
[current_voxel_ID_theta, p1_xy, P_max_pol] = ...
    angular_init(num_polar_sections, sphere_center, sphere_max_radius,...
    r, ray_origin, ray_direction, 'xy');
[current_voxel_ID_phi, p1_xz, P_max_azi] = ...
    angular_init(num_azimuthal_sections, sphere_center, sphere_max_radius,...
    r, ray_origin, ray_direction, 'xz');

% Set initial voxel IDs. 
azimuthal_voxels = [current_voxel_ID_phi];
angular_voxels = [current_voxel_ID_theta];
radial_voxels = [current_voxel_ID_r];

if verbose
    fprintf('\nInitial radial voxel: %d \n', current_voxel_ID_r)
    fprintf('\nInitial theta voxel: %d', current_voxel_ID_theta)
    fprintf('\nInitial phi voxel: %d', current_voxel_ID_phi)
end

% Find the maximum time the ray will be in the grid.
[ta, tb] = radial_intersection_times(ray_origin, ray_direction,...
        sphere_center, sphere_max_radius);
t_grid = max(ta,tb);

% Determine the correct time to begin the traversal phase. If the ray
% starts outside the grid at t_begin, snap to the grid and find the
% corresopnding start time.
if approximatelyEqual(r,sphere_max_radius,1e-12,1e-8)
    if ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
        t_start = (p1_xy(1) - ray_origin(1))/ray_direction(1);
    elseif ~approximatelyEqual(ray_direction(2),0.0,1e-12,1e-8)
        t_start = (p1_xy(2) - ray_origin(2))/ray_direction(2);
    else
        t_start = (p1_xz(2) - ray_origin(3))/ray_direction(3);
    end
else
    t_start = t_begin;
end
    
% TRAVERSAL PHASE
t = t_start;
t_end = min(t_grid, t_end);
time_test = t_end-t_start;
traversal_time = 0;
previous_transition_flag = false;
while t < t_end
    t_past = t;
    % 1. Calculate tMaxR, tMaxTheta, tMaxPhi
    [tMaxR, tStepR, previous_transition_flag] = radial_hit(ray_origin, ...
        ray_direction, current_voxel_ID_r, sphere_center, sphere_max_radius, ...
        delta_radius, t, previous_transition_flag, verbose);
    [tMaxTheta, tStepTheta] = polar_hit(ray_origin, ray_direction, current_voxel_ID_theta,...
        sphere_center, P_max_pol, sphere_max_radius, t, t_end, verbose);
    [tMaxPhi, tStepPhi] = azimuthal_hit(ray_origin, ray_direction, current_voxel_ID_phi,...
        sphere_center, P_max_azi, sphere_max_radius, t, t_end, verbose);   
    % 2. Compare tMaxR, tMaxTheta, tMaxPhi
    if (strictlyLess(tMaxTheta,tMaxR,1e-12,1e-8)  && ...
            strictlyLess(tMaxTheta, tMaxPhi,1e-12,1e-8)) ...
            && strictlyLess(t,tMaxTheta,1e-12,1e-8) && ...
            strictlyLess(tMaxTheta,t_end,1e-12,1e-8) 
        fprintf("theta min")
        % Note 1: Case tMaxTheta is a minimum
        % Note 2: When the ray only intersects one radial shell but crosses an
        % angular boundary, we need the second half of disjunction
        t = tMaxTheta;
        current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_polar_sections);
    elseif strictlyLess(tMaxR, tMaxTheta, 1e-12,1e-8) && ...
            strictlyLess(tMaxR, tMaxPhi, 1e-12,1e-8) && ...
            strictlyLess(t,tMaxR,1e-12,1e-8) && ...
            strictlyLess(tMaxR,t_end, 1e-12,1e-8)
        % Case tMaxR is minimum
        t = tMaxR;
        current_voxel_ID_r = current_voxel_ID_r + tStepR; 
    elseif strictlyLess(tMaxPhi,tMaxTheta,1e-12,1e-8) && ...
            strictlyLess(tMaxPhi,tMaxR,1e-12,1e-8) && ...
            strictlyLess(t,tMaxPhi,1e-12,1e-8) && ...
            strictlyLess(tMaxPhi,t_end,1e-12,1e-8)
        % Case tMaxPhi is minimum
        t = tMaxPhi;
        current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections); 
    elseif approximatelyEqual(tMaxPhi,tMaxTheta,1e-12,1e-8) && ...
            approximatelyEqual(tMaxPhi,tMaxR,1e-12,1e-8) && ...
            strictlyLess(t,tMaxR,1e-12,1e-8) && ...
            strictlyLess(tMaxR,t_end,1e-12,1e-8) 
        % Triple boundary intersection
        t = tMaxPhi;
        current_voxel_ID_r = current_voxel_ID_r + tStepR;
        current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_polar_sections);
        current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);
    elseif approximatelyEqual(tMaxPhi, tMaxTheta, 1e-12,1e-8) && ...
            strictlyLess(t,tMaxPhi,1e-12,1e-8) && ...
            strictlyLess(tMaxPhi,t_end, 1e-12,1e-8)
        % Phi, Theta equal mins
            t = tMaxPhi;
            current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_polar_sections);
            current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);
    elseif approximatelyEqual(tMaxTheta,tMaxR,1e-12,1e-8) && ...
            strictlyLess(t,tMaxR,1e-12,1e-8) && ...
            strictlyLess(tMaxR, t_end,1e-12,1e-8)
        % R, Theta equal mins
            t = tMaxTheta;
            current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_polar_sections);
            current_voxel_ID_r = current_voxel_ID_r + tStepR;
    elseif approximatelyEqual(tMaxR,tMaxPhi,1e-12,1e-8) && ...
            strictlyLess(t,tMaxR,1e-12,1e-8) && ...
            strictlyLess(tMaxR, t_end, 1e-12,1e-8) 
        % R, Phi equal mins
            t = tMaxR;
            current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);
            current_voxel_ID_r = current_voxel_ID_r + tStepR;
    else
       traversal_time = traversal_time + (t_end - t_past);
       return;
    end
    traversal_time = traversal_time + (t - t_past);
    radial_voxels = [radial_voxels, current_voxel_ID_r];
    angular_voxels = [angular_voxels, current_voxel_ID_theta];
    azimuthal_voxels = [azimuthal_voxels, current_voxel_ID_phi];
end

end
