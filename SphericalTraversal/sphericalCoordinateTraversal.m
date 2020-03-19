function [radial_voxels, angular_voxels, phi_voxels] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, sphere_center, ...
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
%    num_azimuthal_sections: The number of azimuthal sections in the sphere.
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
%    phi_voxels: A list of the phi voelx that were hit by the ray.
%    Note: These lists, used in conjunction, will produce the path of the ray
%    through the voxels using each point. For example,
%    [radial_voxels(1), angular_voxels(1)] is the first voxel the ray
%    travels through.
close all;
circle_center_x = circle_center(1);
circle_center_y = circle_center(2);
circle_center_z = circle_center(3);
ray_origin_x = ray_origin(1);
ray_origin_y = ray_origin(2);
ray_origin_z = ray_origin(3);
ray_direction_x = ray_direction(1);
ray_direction_y = ray_direction(2);
ray_direction_z = ray_direction(3);

min_bound_x = min_bound(1);
min_bound_y = min_bound(2);
min_bound_z = min_bound(3);
max_bound_x = max_bound(1);
max_bound_y = max_bound(2);
max_bound_z = max_bound(3);

ray_start = ray_origin + t_begin * ray_direction;
ray_start_x = ray_start(1);
ray_start_y = ray_start(2);
ray_start_z = ray_start(3);

ray_end = ray_origin + t_end * ray_direction;
ray_end_x = ray_end(1);
ray_end_y = ray_end(2);
ray_end_z = ray_end(3);

angular_voxels = [];
radial_voxels = [];
phi_voxels = [];

delta_radius = sphere_max_radius / num_radial_sections;

% INITIALIZATION PHASE
% Determine ray location at t_begin.
p = ray_origin + t_begin.*ray_direction;
ray_sphere_vector = [sphere_center(1) - p(1); sphere_center(2) - p(2); sphere_center(3) - p(3)]';
% Find the radial shell containing the ray at t_begin.
r = delta_radius;
while (ray_sphere_vector(1)^2 + ray_sphere_vector(2)^2 + ray_sphere_vector(3) > r^2) && r < sphere_max_radius
    r = r + delta_radius;
end

% Find the intersection times for the ray and the radial shell containing the ray
% at t_begin; in order to determine if ray intersects grid.
ray_unit_vector = 1 / sqrt(ray_direction(1)^2 + ray_direction(2)^2)...
    .* [ray_direction(1);  ray_direction(2)]';
v = dot(ray_circle_vector,ray_unit_vector);
discr = r^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2);
d = sqrt(discr);
pa = ray_origin + (v-d).*ray_unit_vector;
pb = ray_origin + (v+d).*ray_unit_vector;
tol = 10^-16;
if (abs(ray_direction(1)) < tol)
    t1 = (pa(2) - ray_origin(2))/ray_direction(2);
    t2 = (pb(2) - ray_origin(2))/ray_direction(2);
else
    t1 = (pa(1) - ray_origin(1))/ray_direction(1);
    t2 = (pb(1) - ray_origin(1))/ray_direction(1);
end

% The ray may not intersect the grid at all. 
% In particular, if the ray is outside the grid at t_begin.
if t1 < t_begin && t2 < t_begin 
    if verbose
        fprintf("\nRay does not intersect polar grid for t_begin.")
    end
    return;
end

% It may be a tangent hit
 if abs(t1 - t2) < tol
    if verbose
        fprintf("\nTangent hit.")
    end
    return;
end

% If there is a ray/shell intersection, then set the radial voxel ID. 
current_voxel_ID_r = 1 + (circle_max_radius - r)/delta_radius;
if verbose
    fprintf('RADIAL HIT INITIALIZED.\n')
end

% II. Calculate Voxel ID Theta.
tol = 10^-16;
if abs(ray_origin - circle_center) < tol 
    % If the ray starts at the origin, we need to perturb slightly along its path to find the
    % correct angular voxel
    pert_t = 0.1;
    pert_x = ray_start_x + ray_direction_x * pert_t;
    pert_y = ray_start_y + ray_direction_y * pert_t;
    current_voxel_ID_theta = floor(atan2(pert_y - circle_center_y, pert_x - circle_center_x) * num_angular_sections / (2 * pi));
else
    current_voxel_ID_theta = floor(atan2(ray_start_y - circle_center_y, ray_start_x - circle_center_x) * num_angular_sections / (2 * pi));
end
if current_voxel_ID_theta < 0
    current_voxel_ID_theta = num_angular_sections + current_voxel_ID_theta;
end

angular_voxels = [current_voxel_ID_theta];
radial_voxels = [current_voxel_ID_r];

% Find the maximum time the ray will be in the grid
discr = circle_max_radius^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2);
d = sqrt(discr);
pa_max = ray_origin + (v-d).*ray_unit_vector;
pb_max = ray_origin + (v+d).*ray_unit_vector;
if ray_direction(1) < tol
    t1 = (pa_max(2) - ray_origin(2))/ray_direction(2);
    t2 = (pb_max(2) - ray_origin(2))/ray_direction(2);
else
    t1 = (pa_max(1) - ray_origin(1))/ray_direction(1);
    t2 = (pb_max(1) - ray_origin(1))/ray_direction(1);
end
t_grid = max(t1,t2);

% III. TRAVERSAL PHASE
t = t_begin;
previous_transition_flag = false;
while t < min(t_grid,t_end)
    
    % 1. Calculate tMaxR, tMaxTheta
    [tMaxR, tStepR, previous_transition_flag] = radial_hit(ray_origin, ray_direction, ...
        current_voxel_ID_r, circle_center, circle_max_radius, delta_radius, t, ray_unit_vector, ray_circle_vector, v, previous_transition_flag, verbose);
    [tMaxTheta, tStepTheta] = angular_hit(ray_origin, ray_direction, current_voxel_ID_theta,...
        num_angular_sections, circle_center, t, verbose);
    
    % 2. Compare tMaxTheta, tMaxR
    if (tMaxTheta < tMaxR || current_voxel_ID_r + tStepR == 0) && t < tMaxTheta && tMaxTheta < min(t_grid,t_end)
        % when the ray only intersects one radial shell but crosses an
        % angular boundary, we need the second half of conditional
        t = tMaxTheta;
        current_voxel_ID_theta = current_voxel_ID_theta + tStepTheta;
        if current_voxel_ID_theta < 0
            current_voxel_ID_theta = num_angular_sections + current_voxel_ID_theta;
        end
        if verbose
            new_x_position = ray_origin(1) + ray_direction(1) * tMaxTheta;
            new_y_position = ray_origin(2) + ray_direction(2) * tMaxTheta;
            text(new_x_position, new_y_position, 'POI_t');
            fprintf('ANGULAR HIT.\n')
        end
        angular_voxels = [angular_voxels, current_voxel_ID_theta];
        radial_voxels = [radial_voxels, current_voxel_ID_r];
    elseif tMaxTheta - tMaxR < tol && tMaxR < min(t2,t_end)
        % For the case when the ray simultaneously hits a radial and
        % angular boundary.
        t = tMaxR;
        p = ray_origin + t.*ray_direction;
        current_voxel_ID_r = current_voxel_ID_r + tStepR;
        current_voxel_ID_theta = current_voxel_ID_theta + tStepTheta;
        if current_voxel_ID_theta < 0
            current_voxel_ID_theta = num_angular_sections + current_voxel_ID_theta;
        end
        if verbose
            new_x_position = ray_origin_x + ray_direction_x * tMaxR;
            new_y_position = ray_origin_y + ray_direction_y * tMaxR;       
            if tStepR == 1
                text(new_x_position, new_y_position, 'POI_{rt}');
                fprintf('DOUBLE HIT (inward).\n');
            else
                text(new_x_position, new_y_position, 'POI_{rt}');
                fprintf('DOUBLE HIT (outward).\n');
            end
        end
        angular_voxels = [angular_voxels, current_voxel_ID_theta];
        radial_voxels = [radial_voxels, current_voxel_ID_r];
    elseif tMaxR < min(t2,t_end) && current_voxel_ID_r + tStepR ~= 0
        t = tMaxR;
        p = ray_origin + t.*ray_direction;
        current_voxel_ID_r = current_voxel_ID_r + tStepR;
        if verbose
          new_x_position = ray_origin_x + ray_direction_x * tMaxR;
          new_y_position = ray_origin_y + ray_direction_y * tMaxR;       
          if tStepR == 1
            text(new_x_position, new_y_position, 'POI_r');
            fprintf('RADIAL HIT (inward).\n');
          else
            text(new_x_position, new_y_position, 'POI_r');
            fprintf('RADIAL HIT (outward).\n');
          end
        end
        angular_voxels = [angular_voxels, current_voxel_ID_theta];
        radial_voxels = [radial_voxels, current_voxel_ID_r];
    else
        t = min(t2,t_end);
        if verbose
          fprintf('Grid exit.\n');
        end
    end    
end
end
