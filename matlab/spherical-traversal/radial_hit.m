function[tMaxR, tStepR, transition_flag]=radial_hit(ray_origin, ray_direction, current_radial_voxel, ...
        sphere_center, sphere_max_radius, delta_radius, t, prev_transition_flag, verbose)
% Determines whether a radial hit occurs for the given ray.
% This follows closely the mathematics presented in:
% http://cas.xav.free.fr/Graphics%20Gems%204%20-%20Paul%20S.%20Heckbert.pdf
%
% Input:
%    ray_origin: The origin of the ray.
%    ray_direction: The direction of the ray.
%    current_radial_voxel: The current radial voxel the ray is located in.
%    sphere_center: The center of the circle.
%    sphere_max_radius: The max radius of the circle.
%    delta_radius: The delta of the radial sections.
%    t: The current time parameter of the ray.
%    ray_unit_vector: The ray unit vector.
%    ray_circle_vector: The vector difference between the ray origin and the circle (spherical disc) origin.
%    v: The dot product between the ray_unit_vector and the ray_circle_vector.
%    prev_transition_flag: Determines whether the previous radial traversal
%    was a 'transition.' A transition is defined as the turning point from
%    inward movement to outward movement from the circle's origin. Another
%    way this can be denoted is sequential hits with equal radii.
%    verbose: Determines whether debug print statements are enabled.
%
% Returns:
%    tMaxR: is the time at which a hit occurs for the ray at the next point of intersection.
%    tStepR: The voxel traversal value: 0, +1, -1.
%    transition_flag: Determines whether the current voxel traversal was a
%    'transition.'

if verbose
    fprintf('\n--radial_hit-- \nCurrent Radial Voxel: %d', current_radial_voxel)
end


% Find the "current" radius of the ray based on the current 
% radial voxel.
r = sphere_max_radius - delta_radius * (current_radial_voxel - 1);
% Find the previous radius 
r_a = max(r - delta_radius , delta_radius);
% Find the next radius
% In the case that the ray has sequential hits with equal radii 
% e.g. the innermost radial disc, ensure that proper radii are
% being checked: without this, code skips ahead one radial boundary.
if prev_transition_flag
r_b = min(r, sphere_max_radius);
else
r_b = min(r + delta_radius, sphere_max_radius);
end
% Find the intersection times for the ray and the previous and 
% next radial discs.
% Initialize two arrays for the intersection times of the ray 
% and neighboring (preiovus and next) spherical discs.
time_array_a = [];
[~,~,discra] = radial_intersection_points(ray_origin, ray_direction,...
        sphere_center, r_a);
if (discra >= 0 )        
    [t1,t2] = radial_intersection_times(ray_origin, ray_direction,...
        sphere_center, r_a);
else
    r_a = r_a + delta_radius;
    [t1,t2] = radial_intersection_times(ray_origin, ray_direction,...
        sphere_center, r_a);
end
time_array_a(1) = t1;
time_array_a(2) = t2;

time_array_b = [];
[~,~,discrb] = radial_intersection_points(ray_origin, ray_direction,...
        sphere_center, r_b);
if discrb >= 0      
    [t1,t2] = radial_intersection_times(ray_origin, ray_direction,...
        sphere_center, r_b);
    time_array_b(1) = t1;
    time_array_b(2) = t2;
end
time_array = [time_array_a time_array_b];
time = time_array(time_array > t);

% If it is a "glancing" blow (i.e. ray is tangent to the circle). This
% occurs when the two intersections times are equal.
if length(time)>1 && ...
        approximatelyEqual(time_array(1),time_array(2),1e-12,1e-8)
    tMaxR = time(1);
    p = ray_origin + tMaxR.*ray_direction;
    r_new = sqrt((p(1) - sphere_center(1))^2 + (p(2) - sphere_center(2))^2 + ...
        (p(3) - sphere_center(3))^2);
    tStepR = 0;
    if approximatelyEqual(r,r_new,1e-12,1e-8)
        transition_flag = true;
    else
        transition_flag = false;
    end
    if verbose
        fprintf("\nglancing blow\n")
    end
    return;
end

if isempty(time)
    tMaxR = inf;
    tStepR = 0;
    transition_flag = false;
    if verbose
        fprintf("\nNo intersection for a radial hit for current r: %d", r)
    end
    return;
end

tMaxR = time(1);
p = ray_origin + tMaxR.*ray_direction;
    r_new = sqrt((p(1) - sphere_center(1))^2 + (p(2) - sphere_center(2))^2 + ...
        (p(3) - sphere_center(3))^2);

%Flag for case that the ray has sequential hits with equal radii.
if approximatelyEqual(r,r_new,1e-12,1e-8)
    transition_flag = true;
else
    transition_flag = false;
end

if strictlyLess(r_new,r,1e-12,1e-8)
    tStepR = 1;
else
    tStepR = -1;
end

if verbose
    fprintf('\ntMaxR: %d \n', tMaxR);
end

end
