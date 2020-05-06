function [radial_voxels, angular_voxels, azimuthal_voxels, time_test, traversal_time] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, sphere_center, ...
    sphere_max_radius, num_radial_sections, num_angular_sections, num_azimuthal_sections, t_begin, t_end, verbose)
% Input:
%    min_bound: The lower left corner of the bounding box.
%    max_bound: The upper right corner of the bounding box.
%    ray origin: The origin of the ray in (x, y, z) coordinates.
%    ray direction: The direction of the ray in (x, y, z) coordinates.
%    sphere_center: The x, y location of the center of the sphere.
%    sphere_max_radius: The largest that encompasses the sphere.
%    num_radial_sections: The number of radial sections in the sphere.
%    num_angular_sections: The number of angular sections in the sphere.
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
%    num_angular_sections > 0
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
delta_radius = sphere_max_radius / num_radial_sections;

% Determine ray location at t_begin.
p = ray_origin + t_begin .* ray_direction;
ray_sphere_vector = [sphere_center(1) - p(1); sphere_center(2) - p(2); sphere_center(3) - p(3)]';
% Find the radial shell containing the ray at t_begin.
r = delta_radius;
while (ray_sphere_vector(1)^2 + ray_sphere_vector(2)^2 + ray_sphere_vector(3)^2 > r^2) && r < sphere_max_radius
    r = r + delta_radius;
end

% Find the intersection times for the ray and the radial shell containing 
% the ray at t_begin; in order to determine if the ray intersects the grid.
ray_unit_vector = 1 / sqrt(ray_direction(1)^2 + ray_direction(2)^2 + ray_direction(3)^2)...
    .* [ray_direction(1);  ray_direction(2); ray_direction(3)]';
v = dot(ray_sphere_vector, ray_unit_vector);
discr = r^2 - (dot(ray_sphere_vector,ray_sphere_vector) - v^2);
d = sqrt(discr);
pa = ray_origin + (v-d) .* ray_unit_vector;
pb = ray_origin + (v+d) .* ray_unit_vector;
% Calculate the time of entrance and exit of the ray.
if ~approximatelyEqual(ray_direction(2),0.0,1e-12,1e-8)
    % Use the y-direction if it is non-zero.
    t1 = (pa(2) - ray_origin(2)) / ray_direction(2);
    t2 = (pb(2) - ray_origin(2)) / ray_direction(2);
elseif ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
    % Use the x-direction if it is non-zero.
    t1 = (pa(1) - ray_origin(1)) / ray_direction(1);
    t2 = (pb(1) - ray_origin(1)) / ray_direction(1);
else
    % Use the z-direction if it is non-zero.
    t1 = (pa(3) - ray_origin(3)) / ray_direction(3);
    t2 = (pb(3) - ray_origin(3)) / ray_direction(3);
end

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

% If there is a ray/shell intersection, then set the radial voxel ID. 
current_voxel_ID_r = 1 + (sphere_max_radius - r) / delta_radius;
if verbose
    fprintf('RADIAL HIT INITIALIZED.\n')
end

% Calculate Voxel ID Theta.
delta_theta = 2 * pi/ num_angular_sections;
delta_phi = 2 * pi/ num_azimuthal_sections;
% Create an array of values representing the points of intersection between 
% the lines corresponding to angular voxels boundaries and the initial 
% radial voxel of the ray. Note that spherical coordinates are unnecessary,
% this is just marking the voxel boundaries in the plane. 
i = 1;
k = 0;
trig_ang = zeros(num_angular_sections,2);
while k <= 2*pi 
    trig_ang(i,1) = cos(k);
    trig_ang(i,2) = sin(k);
    k = k + delta_theta;
    i = i + 1;
end
P_max_ang = sphere_max_radius .* trig_ang + [sphere_center(1) sphere_center(2)];
P_ang = r .* trig_ang + [sphere_center(1) sphere_center(2)];

i = 1;
j = 0;
trig_azi = zeros(num_azimuthal_sections,2);
while j <= 2*pi 
    trig_azi(i,1) = cos(j);
    trig_azi(i,2) = sin(j);
    i = i + 1;
    j = j + delta_phi;
end
P_max_azi = sphere_max_radius .* trig_azi + [sphere_center(1) sphere_center(3)];
P_azi = r .* trig_azi + [sphere_center(1) sphere_center(3)];

% Find the point of intersection between the vector created by the
% ray intersection with the initial sphere radius and the sphere center.
ray_origin_xy = [ray_origin(1), ray_origin(2)];
sphere_center_xy = [sphere_center(1), sphere_center(2)];
if approximatelyEqual(ray_origin_xy,sphere_center_xy,1e-12,1e-8)    
    % If the ray starts at the origin, we need to perturb slightly along its path to find the
    % correct angular voxel
    pert_t = 0.1;
    pert_x = ray_origin(1) + ray_direction(1) * pert_t;
    pert_y = ray_origin(2) + ray_direction(2) * pert_t;
    a = sphere_center(1) - pert_x;
    b = sphere_center(2) - pert_y;
elseif approximatelyEqual(r,sphere_max_radius,1e-12,1e-8)
% If the ray origin is outside the grid, this will snap to the grid and use the
% resulting intersection point as the ray origin for this calculation.
    a = sphere_center(1) - pa(1);
    b = sphere_center(2) - pa(2);
else
    a = sphere_center(1) - ray_origin(1);
    b = sphere_center(2) - ray_origin(2);
end
l = sqrt(a^2 + b^2);
if approximatelyEqual(l,0.0,1e-12,1e-8)
    % if l is 0.0 then only the z-dir of the ray is non-zero; set the
    % initial point s.t. voxel_ID_theta init is 0.
    p1_xy = [sphere_center(1) + r, sphere_center(2)];
else
    p1_xy = sphere_center_xy - (r/l) .* [a b];
end
% This point will lie between two angular voxel boundaries iff the angle between
% it and the angular boundary intersection points along the circle of 
% max radius is obtuse. Equality represents the case when the
% point lies on an angular boundary. 
i = 1;    
while i < length(P_ang)
    d1 = (P_ang(i,1)-p1_xy(1))^2 + (P_ang(i,2)-p1_xy(2))^2;
    d2 = (P_ang(i+1,1)-p1_xy(1))^2 + (P_ang(i+1,2)-p1_xy(2))^2;
    d3 = (P_ang(i,1)-P_ang(i+1,1))^2 + (P_ang(i,2)-P_ang(i+1,2))^2;
    if strictlyLess(d1+d2,d3,1e-12,1e-8) || ...
                approximatelyEqual(d1+d2,d3,1e-12,1e-8)
        current_voxel_ID_theta = i - 1;
        i = length(P_ang);
    end
    i = i + 1;
end

ray_origin_xz = [ray_origin(1), ray_origin(3)];
sphere_center_xz = [sphere_center(1), sphere_center(3)];
if approximatelyEqual(ray_origin_xz,sphere_center_xz,1e-12,1e-8)    
    % If the ray starts at the origin, we need to perturb slightly along its path to find the
    % correct angular voxel
    pert_t = 0.1;
    pert_x = ray_origin(1) + ray_direction(1) * pert_t;
    pert_z = ray_origin(3) + ray_direction(3) * pert_t;
    a = sphere_center(1) - pert_x;
    c = sphere_center(3) - pert_z;
elseif approximatelyEqual(r,sphere_max_radius,1e-12,1e-8)
% If the ray origin is outside the grid, this will snap to the grid and use the
% resulting intersection point as the ray origin for this calculation.
    a = sphere_center(1) - pa(1);
    c = sphere_center(3) - pa(3);
else
    a = sphere_center(1) - ray_origin(1);
    c = sphere_center(3) - ray_origin(3);
end
l = sqrt(a^2 + c^2);
if approximatelyEqual(l,0.0,1e-12,1e-8)
    % if l is 0.0 then only the z-dir of the ray is non-zero; set the
    % initial point s.t. voxel_ID_theta init is 0.
    p1_xz = [sphere_center(1) + r, sphere_center(3)];
else
    p1_xz = sphere_center_xz - (r/l) .* [a c];
end
% This point will lie between two angular voxel boundaries iff the angle between
% it and the angular boundary intersection points along the circle of 
% max radius is obtuse. Equality represents the case when the
% point lies on an angular boundary. 
i = 1;    
while i < length(P_azi)
    d1 = (P_azi(i,1)-p1_xz(1))^2 + (P_azi(i,2)-p1_xz(2))^2;
    d2 = (P_azi(i+1,1)-p1_xz(1))^2 + (P_azi(i+1,2)-p1_xz(2))^2;
    d3 = (P_azi(i,1)-P_azi(i+1,1))^2 + (P_azi(i,2)-P_azi(i+1,2))^2;
    if strictlyLess(d1+d2,d3,1e-12,1e-8) || ...
                approximatelyEqual(d1+d2,d3,1e-12,1e-8)
        current_voxel_ID_phi = i - 1;
        i = length(P_azi);
    end
    i = i + 1;
end

azimuthal_voxels = [current_voxel_ID_phi];
angular_voxels = [current_voxel_ID_theta];
radial_voxels = [current_voxel_ID_r];

if verbose
    fprintf('\nInitial phi voxel: %d', current_voxel_ID_phi)
    fprintf('\nInitial theta voxel: %d', current_voxel_ID_theta)
    fprintf('\nInitial radial voxel: %d \n', current_voxel_ID_r)
end

% Find the maximum time the ray will be in the grid.
discr = sphere_max_radius^2 - (dot(ray_sphere_vector,ray_sphere_vector) - v^2);
d = sqrt(discr);
pa_max = ray_origin + (v-d) .* ray_unit_vector;
pb_max = ray_origin + (v+d) .* ray_unit_vector;
if ~approximatelyEqual(ray_direction(2),0.0,1e-12,1e-8)
    t1 = (pa_max(2) - ray_origin(2)) / ray_direction(2);
    t2 = (pb_max(2) - ray_origin(2)) / ray_direction(2);
elseif ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
    t1 = (pa_max(1) - ray_origin(1)) / ray_direction(1);
    t2 = (pb_max(1) - ray_origin(1)) / ray_direction(1);
else 
    t1 = (pa_max(3) - ray_origin(3)) / ray_direction(3);
    t2 = (pb_max(3) - ray_origin(3)) / ray_direction(3);
end
t_grid = max(t1,t2);

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
    [tMaxR, tStepR, previous_transition_flag] = radial_hit(ray_origin, ray_direction, ...
        current_voxel_ID_r, sphere_center, sphere_max_radius, delta_radius, t, ray_unit_vector, ...
        ray_sphere_vector, v, previous_transition_flag, verbose);
    [tMaxTheta, tStepTheta] = angular_hit(ray_origin, ray_direction, current_voxel_ID_theta,...
        sphere_center, P_max_ang, sphere_max_radius, t, t_end, verbose);
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
        current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_angular_sections);
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
        current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_angular_sections);
        current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);
    elseif approximatelyEqual(tMaxPhi, tMaxTheta, 1e-12,1e-8) && ...
            strictlyLess(t,tMaxPhi,1e-12,1e-8) && ...
            strictlyLess(tMaxPhi,t_end, 1e-12,1e-8)
        % Phi, Theta equal mins
            t = tMaxPhi;
            current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_angular_sections);
            current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);
    elseif approximatelyEqual(tMaxTheta,tMaxR,1e-12,1e-8) && ...
            strictlyLess(t,tMaxR,1e-12,1e-8) && ...
            strictlyLess(tMaxR, t_end,1e-12,1e-8)
        % R, Theta equal mins
            t = tMaxTheta;
            current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_angular_sections);
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
