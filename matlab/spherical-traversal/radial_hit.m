function[tMaxR, tStepR, transition_flag]=radial_hit(ray_origin, ray_direction, current_radial_voxel, ...
        sphere_center, sphere_max_radius, delta_radius, t, ray_unit_vector, ray_circle_vector, v, prev_transition_flag, verbose)
% Determines whether a radial hit occurs for the given ray.
% This follows closely the mathematics presented in:
% http://cas.xav.free.fr/Graphics%20Gems%204%20-%20Paul%20S.%20Heckbert.pdf
%
% Note also we've currently provided a tolerance 'tol' below. This is used
% to account for floating point error.
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
time_array_b = [];
tol = 10^-13;
discr = r_a^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2);
if (discr >= 0 )        
    d = sqrt(discr);
    ta = (v-d);
    tb = (v+d);
    pa = ray_origin + ta.*ray_unit_vector;
    pb = ray_origin + tb.*ray_unit_vector;    
    if abs(ray_direction(2)) > tol
        t1 = (pa(2) - ray_origin(2))/ray_direction(2);
        t2 = (pb(2) - ray_origin(2))/ray_direction(2);
    elseif abs(ray_direction(1)) > tol
        t1 = (pa(1) - ray_origin(1))/ray_direction(1);
        t2 = (pb(1) - ray_origin(1))/ray_direction(1);
    else
        t1 = (pa(3) - ray_origin(3))/ray_direction(3);
        t2 = (pb(3) - ray_origin(3))/ray_direction(3);
    end
    time_array_a(1) = t1;
    time_array_a(2) = t2;
else
    r_a = r_a + delta_radius;
    discr = r_a^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2);
    d = sqrt(discr);
    ta = (v-d);
    tb = (v+d);
    pa = ray_origin + ta.*ray_unit_vector;
    pb = ray_origin + tb.*ray_unit_vector;
    if abs(ray_direction(2)) > tol
        t1 = (pa(2) - ray_origin(2))/ray_direction(2);
        t2 = (pb(2) - ray_origin(2))/ray_direction(2);
    elseif abs(ray_direction(1)) > tol
        t1 = (pa(1) - ray_origin(1))/ray_direction(1);
        t2 = (pb(1) - ray_origin(1))/ray_direction(1);
    else
        t1 = (pa(3) - ray_origin(3))/ray_direction(3);
        t2 = (pb(3) - ray_origin(3))/ray_direction(3);
    end
    time_array_a(1) = t1;
    time_array_a(2) = t2;
end
discr = r_b^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2);
if discr >= 0      
    d = sqrt(discr);
    ta = (v-d);
    tb = (v+d);
    pa = ray_origin + ta.*ray_unit_vector;
    pb = ray_origin + tb.*ray_unit_vector;
    if abs(ray_direction(2)) > tol
        t1 = (pa(2) - ray_origin(2))/ray_direction(2);
        t2 = (pb(2) - ray_origin(2))/ray_direction(2);
    elseif abs(ray_direction(1)) > tol
        t1 = (pa(1) - ray_origin(1))/ray_direction(1);
        t2 = (pb(1) - ray_origin(1))/ray_direction(1);
    else
        t1 = (pa(3) - ray_origin(3))/ray_direction(3);
        t2 = (pb(3) - ray_origin(3))/ray_direction(3);
    end
    time_array_b(1) = t1;
    time_array_b(2) = t2;
end
time_array = [time_array_a time_array_b];
time = time_array(time_array > t);

% If it is a "glancing" blow (i.e. ray is tangent to the circle). This
% occurs when the two intersections times are equal.
if length(time)>1 &&  abs(time_array(1)-time_array(2)) < tol
    tMaxR = time(1);
    p = ray_origin + tMaxR.*ray_direction;
    r_new = sqrt((p(1) - sphere_center(1))^2 + (p(2) - sphere_center(2))^2 + ...
        (p(3) - sphere_center(3))^2);
    tStepR = 0;
    if (abs(r - r_new) < tol)
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
if (abs(r - r_new) < tol)
    transition_flag = true;
else
    transition_flag = false;
end

if r_new - r < 0 && abs(r_new - r) > tol && ~(abs(r_new-r)) < tol
    tStepR = 1;
else
    tStepR = -1;
end

if verbose
    fprintf('\ntMaxR: %d \n', tMaxR);
end

end
