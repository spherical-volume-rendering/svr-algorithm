function[tMaxR]=radial_hit_dotprod(ray_origin, ray_direction, current_radial_voxel, ...
        circle_center, circle_max_radius, delta_radius, t, verbose)
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
r = circle_max_radius - delta_radius*(current_radial_voxel - 1);

% Check for intersections with relevant radial neighbors.
ray_unit_vector = 1/sqrt(ray_direction(1)^2 + ray_direction(2)^2)...
    .* [ray_direction(1);  ray_direction(2)]';
ray_circle_vector = [circle_center(1) - ray_origin(1); circle_center(2) - ray_origin(2)]';
v = dot(ray_circle_vector,ray_unit_vector);
r_a = max(r - delta_radius , delta_radius);
r_b = min(r + delta_radius, circle_max_radius);
if (abs(r_a - r_b) < eps)
    r_b = r + 2 * delta_radius;
end

time_array_a = [];
time_array_b = [];
discr = r_a^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2);
if (discr >= 0 )        
    d = sqrt(discr);
    ta = (v-d);
    tb = (v+d);
    pa = ray_origin + ta.*ray_unit_vector;
    pb = ray_origin + tb.*ray_unit_vector;
    t1 = (pa(1) - ray_origin(1))/ray_direction(1);
    t2 = (pb(1) - ray_origin(1))/ray_direction(1);
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
    t1 = (pa(1) - ray_origin(1))/ray_direction(1);
    t2 = (pb(1) - ray_origin(1))/ray_direction(1);
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
    t1 = (pa(1) - ray_origin(1))/ray_direction(1);
    t2 = (pb(1) - ray_origin(1))/ray_direction(1);
    time_array_b(1) = t1;
    time_array_b(2) = t2;
end
time_array = [time_array_a time_array_b];
time = time_array(time_array > t);

if (isempty(time))
    tMaxR = inf;
    if verbose
        fprintf("\nNo intersection for a radial hit for current r: %d", r);
    end
    return;
end
tMaxR = time(1);

if verbose
    fprintf('\ntMaxR: %d \n', tMaxR);
end
end
