function[tMaxR, tStepR, transitionFlag]=radial_hit(ray_origin, ray_direction, current_radial_voxel, ...
        circle_center, circle_max_radius, delta_radius, t, prev_transitionFlag, verbose)
% Determines whether a radial hit occurs for the given ray.
% Input:
%    ray_origin: The origin of the ray.
%    ray_direction: The direction of the ray.
%    current_radial_voxel: The current radial voxel the ray is located in.
%    circle_center: The center of the circle.
%    circle_max_radius: The max radius of the circle.
%    delta_radius: The delta of the radial sections.
%    t: The current time parameter of the ray.
%    verbose: Determines whether debug print statements are enabled.
%
% Returns:
%    tMaxR: is the time at which a hit occurs for the ray at the next point of intersection.

if verbose
    fprintf('\n--radial_hit-- \nCurrent Radial Voxel: %d', current_radial_voxel)
end


% Recapture the current radius of ray.
r = circle_max_radius - delta_radius * (current_radial_voxel - 1)

% Check for intersections with relevant radial neighbors.

% TODO: ray_unit_vector, ray_circle_vector, v can be calculated during
% initialization phase.
ray_unit_vector = 1/sqrt(ray_direction(1)^2 + ray_direction(2)^2)...
    .* [ray_direction(1);  ray_direction(2)]';
ray_circle_vector = [circle_center(1) - ray_origin(1); circle_center(2) - ray_origin(2)]';
v = dot(ray_circle_vector,ray_unit_vector);
r_a = max(r - delta_radius , delta_radius);
% In the case that the ray has sequential hits with equal radii ensure that proper radii are
% being checked: without this, code skips ahead one radial boundary.
if prev_transitionFlag
r_b = min(r, circle_max_radius);
else
r_b = min(r + delta_radius, circle_max_radius);
end

tol = 10e-15;

time_array_a = [];
time_array_b = [];
discr = r_a^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2)
if (discr >= 0 )        
    d = sqrt(discr);
    ta = (v-d);
    tb = (v+d);
    pa = ray_origin + ta.*ray_unit_vector
    pb = ray_origin + tb.*ray_unit_vector
    
    if (ray_direction(1) < tol)
        t1 = (pa(2) - ray_origin(2))/ray_direction(2);
        t2 = (pb(2) - ray_origin(2))/ray_direction(2);
    else
        t1 = (pa(1) - ray_origin(1))/ray_direction(1);
        t2 = (pb(1) - ray_origin(1))/ray_direction(1);
    end
    time_array_a(1) = t1;
    time_array_a(2) = t2;
    
else
    r_a = r_a + delta_radius
    discr = r_a^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2);
    d = sqrt(discr);
    ta = (v-d);
    tb = (v+d);
    pa = ray_origin + ta.*ray_unit_vector;
    pb = ray_origin + tb.*ray_unit_vector;
    if (ray_direction(1) < tol)
        t1 = (pa(2) - ray_origin(2))/ray_direction(2);
        t2 = (pb(2) - ray_origin(2))/ray_direction(2);
    else
        t1 = (pa(1) - ray_origin(1))/ray_direction(1);
        t2 = (pb(1) - ray_origin(1))/ray_direction(1);
    end
    time_array_a(1) = t1;
    time_array_a(2) = t2;
end

discr = r_b^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2);
if (discr >= 0 )        
    d = sqrt(discr);
    ta = (v-d);
    tb = (v+d);
    pa = ray_origin + ta.*ray_unit_vector;
    pb = ray_origin + tb.*ray_unit_vector;
    if (ray_direction(1) < tol)
        t1 = (pa(2) - ray_origin(2))/ray_direction(2);
        t2 = (pb(2) - ray_origin(2))/ray_direction(2);
    else
        t1 = (pa(1) - ray_origin(1))/ray_direction(1);
        t2 = (pb(1) - ray_origin(1))/ray_direction(1);
    end
    time_array_b(1) = t1;
    time_array_b(2) = t2;
end
time_array = [time_array_a time_array_b]
time = time_array(time_array > t)

if (isempty(time))
    tMaxR = inf;
    tStepR = 0;
    transitionFlag = false;
    if verbose
        fprintf("\nNo intersection for a radial hit for current r: %d", r);
    end
    return;
end

%If it is a "glancing" blow (i.e. ray is tangent to the circle)
if (length(time)>1 && abs(time(1)-time(2)) < tol)
    tMaxR = time(1);
    p = ray_origin + tMaxR.*ray_direction;
    r_new = sqrt((p(1) - circle_center(1))^2 + (p(2) - circle_center(2))^2)
    tStepR = 0;
    if (abs(r - r_new) < tol)
        transitionFlag = true;
    else
        transitionFlag = false;
    end
    return;
end

tMaxR = time(1);
p = ray_origin + tMaxR.*ray_direction;
r_new = sqrt((p(1) - circle_center(1))^2 + (p(2) - circle_center(2))^2)

%Flag for case that the ray has sequential hits with equal radii.
if (abs(r - r_new) < tol)
    transitionFlag = true;
else
    transitionFlag = false;
end

if (r_new - r < 0 && abs(r_new - r) > tol && ~(abs(r_new-r)) < tol) 
    tStepR = 1;
else
    tStepR = -1;
end

if verbose
    fprintf('\ntMaxR: %d \n', tMaxR);
end

end
