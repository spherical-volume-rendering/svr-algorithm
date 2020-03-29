function [tMaxR, tStepR] = radial_hit_jenkyR(ray_origin, ray_direction, ...
current_radial_voxel, circle_center, ...
circle_max_radius, delta_radius, jenkyR, t, verbose)
% Determines whether a radial hit occurs for the given ray.
% Input:
%    ray_origin: The origin of the ray.
%    ray_direction: The direction of the ray.
%    current_radial_voxel: The current radial voxel the ray is located in.
%    circle_center: The center of the circle.
%    circle_max_radius: The max radius of the circle.
%    jenkyR: TODO.
%    delta_radius: The delta of the radial sections.
%    t: The current time parameter of the ray.
%    verbose: Determines whether debug print statements are enabled.
%
% Returns:
%    tMaxR: is the time at which a hit occurs for the ray at the next point of intersection.
%    tStepR: The direction of step into the next radial voxel: 0, +1, -1.
ray_direction_x = ray_direction(1);
ray_direction_y = ray_direction(2);
circle_center_x = circle_center(1);
circle_center_y = circle_center(2);
ray_origin_x = ray_origin(1);
ray_origin_y = ray_origin(2);

if verbose
    fprintf('\n--radial_hit-- \nCurrent Radial Voxel: %d\n', current_radial_voxel);
end

% (1)   (x - circle_center_x)^2 + (y - circle_center_y)^2 = current_radius^2
% (2)    x = ray_origin_x + ray_direction_x(t)
% (3)    y = ray_origin_y + ray_direction_y(t)
% Plug in x, y in equation (1), then solve for t.
% To get point of intersection, plug t back in parametric equation of a ray.
syms cT; % current time
current_radius = circle_max_radius - (abs(delta_radius) * (current_radial_voxel - 1));
cr = (current_radius - jenkyR)^2;

if (cr == 0), cr = current_radius^2; end

intersections_t = solve( ...
    (ray_origin_x + ray_direction_x * cT - circle_center_x)^2 + ...
    (ray_origin_y + ray_direction_y * cT - circle_center_y)^2 ...
    - cr == 0, cT);
sect = double(subs(intersections_t));

if isempty(sect) || imag(sect(1)) ~= 0 || imag(sect(2)) ~= 0
    % Case where no intersection between ray and new radial voxel
    % occurs.
    tStepR = -inf;
    tMaxR = -inf;
    if verbose
        fprintf("No intersection occurs.\n");
    end
    return;
end

if (t <  min(double(subs(intersections_t))))
    tMaxR = min(double(subs(intersections_t)));
    tStepR = 1;
else
    tMaxR = max(double(subs(intersections_t)));
    tStepR = -1;
end
end
