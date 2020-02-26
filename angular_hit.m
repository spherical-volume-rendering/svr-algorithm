function [tMaxTheta, tStepTheta] = angular_hit(ray_origin, ray_direction, current_voxel_ID_theta,...
    num_angular_sections, circle_center, t, verbose)
% Determines whether an angular hit occurs for the given ray.
% Input:
%    ray_origin: vector of the origin of the ray in cartesian coordinate
%    ray_direction: vector of the direction of the ray in cartesian
%                   coordinate
%    current_voxel_ID_theta: the (angular) ID of current voxel
%    num_angular_sections: number of total angular sections on the grid
% Returns:
%    is_angular_hit: true if an angular crossing has occurred, false otherwise.
%    tMaxTheta: is the time at which a hit occurs for the ray at the next point of intersection.
%    tStepTheta: The direction the theta voxel steps. +1, -1, or 0.
if verbose
    fprintf("\n-- angular_hit --")
end

% First calculate the angular interval that current voxID corresponds
% to
delta_theta = 2 * pi / num_angular_sections;
interval_theta = [current_voxel_ID_theta * delta_theta, (current_voxel_ID_theta + 1) * delta_theta];

if verbose
    fprintf("\nCurrent Voxel ID Theta: %d\n", current_voxel_ID_theta)
end


% Calculate the x and y components that correspond to the angular
% boundary for the angular interval
xmin = cos(min(interval_theta));
xmax = cos(max(interval_theta));
ymin = sin(min(interval_theta));
ymax = sin(max(interval_theta));

if (single(tan(min(interval_theta))) == ray_direction(2)/ray_direction(1) ...
        && single(tan(max(interval_theta))) == ray_direction(2)/ray_direction(1))
    fprintf("parallel")
    tMaxTheta = inf;
    tStepTheta = 0;
    return;
end

% Solve the systems Az=b to check for intersection
Amin = [xmin, -ray_direction(1); ymin, -ray_direction(2)];
Amax = [xmax, -ray_direction(1); ymax, -ray_direction(2)];
b = [ray_origin(1)-circle_center(1), ray_origin(2)-circle_center(2)]';
if (single(tan(min(interval_theta))) == ray_direction(2)/ray_direction(1))
    zmax = Amax\b;
    zmin = [0 ; 0];
elseif (single(tan(max(interval_theta))) == ray_direction(2)/ray_direction(1))
    zmin = Amin\b;
    zmax = [0 ; 0];
else
    zmin = Amin\b;
    zmax = Amax\b;
end

% We need the radius (r = z[1]) and time (t = z[2]) to be positive or
% else the intersection is null. 
% We need the time step of the traversal to be less than t = zmax(2) or 
% else the intersection is null. 

if (((zmin(1) < 0 || zmin(2) < 0) && (zmax(1) < 0 || zmax(2) < 0)) || (t >= max(zmax(2),zmin(2))))
    tMaxTheta = inf;
    tStepTheta = 0;
    if verbose
        fprintf("angular intersection is null")
    end
    return;
end

% If we hit the min boundary then we decrement theta, else increment;
% assign tMaxTheta
if zmin(1) > 0 && zmin(2)>t
    tStepTheta = -1;
    tMaxTheta = zmin(2);
    if verbose
        fprintf("Hit min angular bound\n");
    end
else
    tStepTheta = 1;
    tMaxTheta = zmax(2);
        if verbose
        fprintf("Hit max angular bound\n");
    end
end

if verbose
    fprintf(['\ntMaxTheta: %d \n' ...
        'tStepTheta: %d \n'], tMaxTheta, tStepTheta);
end
end
