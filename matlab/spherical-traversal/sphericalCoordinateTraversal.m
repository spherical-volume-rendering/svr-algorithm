function [radial_voxels, angular_voxels, azimuthal_voxels] = sphericalCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, sphere_center, ...
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
%
%    Note: These lists, used in conjunction, will produce the path of the ray
%    through the voxels using each voxel transition. For example,
%    [radial_voxels(1), angular_voxels(1), azimuthal_voxels(1)]
%    is the first voxel the ray travels through. If the next traversal is a radial hit,
%    angular_voxels(2) and azimuthal_voxels(2) will remain the same.
close all;
sphere_center_x = sphere_center(1);
sphere_center_y = sphere_center(2);
sphere_center_z = sphere_center(3);
ray_direction_x = ray_direction(1);
ray_direction_y = ray_direction(2);
ray_direction_z = ray_direction(3);

ray_start = ray_origin + t_begin * ray_direction;
ray_start_x = ray_start(1);
ray_start_y = ray_start(2);
ray_start_z = ray_start(3);

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
tol = 10^-16;
% Calculate the time of entrance and exit of the ray.
if (abs(ray_direction(2)) > tol)
    % Use the y-direction if it is non-zero.
    t1 = (pa(2) - ray_origin(2)) / ray_direction(2);
    t2 = (pb(2) - ray_origin(2)) / ray_direction(2);
elseif (abs(ray_direction(1)) > tol)
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
    if verbose
        fprintf("\nRay does not intersect spherical grid for t_begin.")
    end
    return;
end

% It may be a tangent hit.
 if abs(t1 - t2) < tol
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
tol = 10^-16; % Tolerance, to account for floating point rounding error.
if abs(ray_origin - sphere_center) < tol 
    % If the ray starts at the origin, we need to perturb slightly along 
    % its path to find the correct angular voxel.
    pert_t = 0.1;
    pert_x = ray_start_x + ray_direction_x * pert_t;
    pert_y = ray_start_y + ray_direction_y * pert_t;
    current_voxel_ID_theta = floor(atan2(pert_y - sphere_center_y, pert_x - sphere_center_x) * num_angular_sections / (2 * pi));
else
    current_voxel_ID_theta = floor(atan2(ray_start_y - sphere_center_y, ray_start_x - sphere_center_x) * num_angular_sections / (2 * pi));
end
if current_voxel_ID_theta < 0
    current_voxel_ID_theta = num_angular_sections + current_voxel_ID_theta;
end

% Calculate Voxel ID phi.
if abs(ray_origin - sphere_center) < tol 
    % If the ray starts at the origin, we need to perturb slightly along
    % its path to find the correct azimuthal voxel.
    pert_t = 0.1;
    pert_x = ray_start_x + ray_direction_x * pert_t;
    pert_z = ray_start_z + ray_direction_z * pert_t;
    current_voxel_ID_phi = floor(atan2(pert_z - sphere_center_z, pert_x - sphere_center_x) * num_azimuthal_sections / (2 * pi));
else
    current_voxel_ID_phi = floor(atan2(ray_start_z - sphere_center_z, ray_start_x - sphere_center_x) * num_azimuthal_sections / (2 * pi));
end
if current_voxel_ID_phi < 0
    current_voxel_ID_phi = num_azimuthal_sections + current_voxel_ID_phi;
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
if ray_direction(2) > tol
    t1 = (pa_max(2) - ray_origin(2)) / ray_direction(2);
    t2 = (pb_max(2) - ray_origin(2)) / ray_direction(2);
elseif ray_direction(1) > tol
    t1 = (pa_max(1) - ray_origin(1)) / ray_direction(1);
    t2 = (pb_max(1) - ray_origin(1)) / ray_direction(1);
else 
    t1 = (pa_max(3) - ray_origin(3)) / ray_direction(3);
    t2 = (pa_max(3) - ray_origin(3)) / ray_direction(3);
end
t_grid = max(t1,t2);

% TRAVERSAL PHASE
t = t_begin;
t_end = min(t_grid, t_end);
previous_transition_flag = false;

while t < t_end    
    % 1. Calculate tMaxR, tMaxTheta, tMaxPhi
    [tMaxR, tStepR, previous_transition_flag] = radial_hit(ray_origin, ray_direction, ...
        current_voxel_ID_r, sphere_center, sphere_max_radius, delta_radius, t, ray_unit_vector, ...
        ray_sphere_vector, v, previous_transition_flag, verbose);
    [tMaxTheta, tStepTheta] = angular_hit(ray_origin, ray_direction, current_voxel_ID_theta,...
        num_angular_sections, sphere_center, t, verbose);
    [tMaxPhi, tStepPhi] = azimuthal_hit(ray_origin, ray_direction, current_voxel_ID_phi,...
      num_azimuthal_sections, sphere_center, t, verbose);
    
    rStepViolation = current_voxel_ID_r + tStepR == 0;
  
    % 2. Compare tMaxR, tMaxTheta, tMaxPhi
    if ((tMaxTheta < tMaxR && tMaxR < tMaxPhi) || rStepViolation) ...
            && t < tMaxTheta && tMaxTheta < t_end
        % Note 1: Case tMaxTheta is a minimum
        % Note 2: When the ray only intersects one radial shell but crosses an
        % angular boundary, we need the second half of disjunction
        t = tMaxTheta;
        current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_angular_sections);
    elseif tMaxR < tMaxTheta && tMaxR < tMaxPhi && t < tMaxR && tMaxR < t_end ...
            && ~rStepViolation
        % Case tMaxR is minimum
        t = tMaxR;
        current_voxel_ID_r = current_voxel_ID_r + tStepR;   
    elseif tMaxPhi < tMaxTheta && tMaxPhi < tMaxR && t < tMaxPhi && tMaxPhi < t_end
        % Case tMaxPhi is minimum
        t = tMaxPhi;
        current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);   
    elseif abs(tMaxPhi - tMaxTheta) < tol && abs(tMaxPhi - tMaxR) < tol && ...
            t < tMaxR && tMaxR < t_end && ~rStepViolation
        % Triple boundary intersection
        t = tMaxPhi;
        current_voxel_ID_r = current_voxel_ID_r + tStepR;
        current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_angular_sections);
        current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);   
    elseif abs(tMaxPhi - tMaxTheta) < tol && t < tMaxPhi && tMaxPhi < t_end
        % Phi, Theta equal
        if tMaxR < tMaxPhi && t < tMaxR && ~rStepViolation
            % R min
            t = tMaxR;
            current_voxel_ID_r = current_voxel_ID_r + tStepR;
        else
            t = tMaxPhi;
            current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_angular_sections);
            current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);
        end
    elseif abs(tMaxTheta - tMaxR) < tol && t < tMaxR && tMaxR < t_end && ...
            ~rStepViolation
        % R, Theta equal
        if tMaxPhi < tMaxTheta && t < tMaxPhi
            % Phi min
            t = tMaxPhi;
            current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);
        else
            t = tMaxTheta;
            current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_angular_sections);
            current_voxel_ID_r = current_voxel_ID_r + tStepR;
        end
    elseif abs(tMaxR - tMaxPhi) < tol && t < tMaxR && tMaxR < t_end && ...
            ~rStepViolation
        % R, Phi equal
        if tMaxTheta < tMaxR && t < tMaxTheta
            % Theta min
            t = tMaxTheta;
            current_voxel_ID_theta = mod(current_voxel_ID_theta + tStepTheta, num_angular_sections);
        else
            t = tMaxR;
            current_voxel_ID_phi = mod(current_voxel_ID_phi + tStepPhi, num_azimuthal_sections);
            current_voxel_ID_r = current_voxel_ID_r + tStepR;
        end
    else
       return;
    end    
    radial_voxels = [radial_voxels, current_voxel_ID_r];
    angular_voxels = [angular_voxels, current_voxel_ID_theta];
    azimuthal_voxels = [azimuthal_voxels, current_voxel_ID_phi];
end
end
