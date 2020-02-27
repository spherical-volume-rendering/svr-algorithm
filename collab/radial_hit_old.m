function[tMaxR, tStepR, is_radial_hit]=radial_hit_old(ray_origin, ray_direction, ...
    current_radial_voxel, circle_center, ...
    circle_max_radius, delta_radius, r, t, verbose)
%[%is_radial_hit,
... tMaxR, %tStepR
    ... ]
    
% Determines whether a radial hit occurs for the given ray.
% Input:
%    ray_origin: The origin of the ray.
%    ray_direction: The direction of the ray.
%    current_radial_voxel: The current radial voxel the ray is located in.
%    circle_center: The center of the circle.
%    circle_max_radius: The max radius of the circle.
%    num_radial_sections: The number of radial sections.
%    delta_radius: The delta of the radial sections.
%
% Returns:
%    is_radial_hit: true if a radial crossing has occurred, false otherwise.
%    tMaxR: is the time at which a hit occurs for the ray at the next point of intersection.
%    tStepR: The direction of step into the next radial voxel, 0, +1, -1.
ray_direction_x = ray_direction(1);
ray_direction_y = ray_direction(2);
circle_center_x = circle_center(1);
circle_center_y = circle_center(2);
ray_origin_x = ray_origin(1);
ray_origin_y = ray_origin(2);

if verbose
    fprintf('\n--radial_hit-- \nCurrent Radial Voxel: %d', current_radial_voxel)
end

% (1)   (x - circle_center_x)^2 + (y - circle_center_y)^2 = current_radius^2
% (2)    x = ray_origin_x + ray_direction_x(t)
% (3)    y = ray_origin_y + ray_direction_y(t)
% Plug in x, y in equation (1), then solve for t.
% To get point of intersection, plug t back in parametric equation of a ray.
p0 = ray_origin + t.*ray_direction
r = max(r, delta_radius)
ray_unit_vector = 1/sqrt(ray_direction(1)^2 + ray_direction(2)^2)...
    .* [ray_direction(1);  ray_direction(2)]';
ray_circle_vector = [circle_center(1) - ray_origin(1); circle_center(2) - ray_origin(2)]';
v = dot(ray_circle_vector,ray_unit_vector);
discr = r^2 - (dot(ray_circle_vector,ray_circle_vector) - v^2)
if (discr < 0 ) 
    tMaxR = inf;
    tStepR = 0;
    return
end
d = sqrt(discr)
if ((v - d) >= 0 && (v + d) >= 0)
    ta = (v-d)
    tb = (v+d)
    pa = ray_origin + ta.*ray_unit_vector
    pb = ray_origin + tb.*ray_unit_vector
    t1= (pa(1) - ray_origin(1))/ray_direction(1)
    t2 = (pb(1) - ray_origin(1))/ray_direction(1)
    if (t1 > t && t2 > t)
    	tMaxR = t1
        p = pa
    elseif (~(t1 > t) && (t2 > t))
        tMaxR = t2
        p = pb
    end
    distance_to_center1 = (circle_center(1)-p0(1))^2 + (circle_center(2)-p0(2))^2
    distance_to_center2 = (circle_center(1)-p(1))^2 + (circle_center(2)-p(2))^2
    distance_to_center1 - distance_to_center2
    tol = 10^-6
    if (distance_to_center1 - distance_to_center2 > tol) 
       tStepR = 1;
    elseif (distance_to_center2 - distance_to_center1 > tol)
              tStepR = -1;
    else
        tStepR = 0;
    end
else
    tMaxR = inf;
    tStepR = 0;
end


if verbose
    fprintf(['\ntMaxR: %d \n' ...
        'tStepR: %d \n'], tMaxR, tStepR);
end
end
