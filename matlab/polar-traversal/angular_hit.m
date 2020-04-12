function [tMaxTheta, tStepTheta] = angular_hit(ray_origin, ray_direction, current_voxel_ID_theta,...
    num_angular_sections, circle_center, pointArray, t, t_end, verbose)
% Determines whether an angular hit occurs for the given ray.
% Input:
%    ray_origin: vector of the origin of the ray in cartesian coordinate.
%    ray_direction: vector of the direction of the ray in cartesian coordinate.
%    current_voxel_ID_theta: the (angular) ID of current voxel.
%    num_angular_sections: number of total angular sections on the grid.
%    circle_center: The center of the circle.
%    pointArray:
%    t: The current time parameter of the ray.
%    t_end:
%    verbose: Determines whether debug print statements are enabled.
%
% Returns:
%    tMaxTheta: is the time at which a hit occurs for the ray at the next point of intersection.
%    tStepTheta: The direction the theta voxel steps. +1, -1, or 0.
if verbose
    fprintf("\n-- angular_hit --")
end

if verbose
    fprintf("\nCurrent Voxel ID Theta: %d\n", current_voxel_ID_theta)
end

% Get current ray location to use as start point
p = ray_origin + ray_direction * t;
% Get ray exit point to use as end point
p_f = ray_origin + ray_direction * t_end;
% Caculate the ray segment vector
v = p_f - p;
% Calculate the voxel boundary vectors
u_min = circle_center - pointArray(1,:);
u_max = circle_center - pointArray(2,:);
% http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
w_min = pointArray(1,:) - p;
w_max = pointArray(2,:) - p;
perpuv_min = u_min(1) * v(2) - u_min(2) * v(1);
perpuv_max = u_max(1) * v(2) - u_max(2) * v(1);
perpuw_min = u_min(1) * w_min(2) - u_min(2) * w_min(1);
perpuw_max = u_max(1) * w_max(2) - u_max(2) * w_max(1);
perpvw_min = v(1) * w_min(2) - v(2) * w_min(1);
perpvw_max = v(1) * w_max(2) - v(2) * w_max(1);

% check to see if ray is parallel to each angular voxel boundary
% if the ray is parallel to or collinear with a given voxel boundary,
% then intersection does not occur between the given voxel boundary and
% the ray.
if approximatelyEqual(perpuv_min,0.0,1e-12,1e-8) ; parallel_min = 1; else; parallel_min = 0; end
if approximatelyEqual(perpuv_max,0.0,1e-12,1e-8) ; parallel_max = 1; else; parallel_max = 0; end

if parallel_min == 0
    a = perpvw_min / perpuv_min;
    b = perpuw_min / perpuv_min;
    if (a < 0 || a > 1) || (b < 0 || b > 1)
        intersect_min = 0;
    else
        intersect_min = 1;
        p_min = p  + b * v;
    end
else
    intersect_min = 0;
end

if parallel_max == 0
    a = perpvw_max / perpuv_max;
    b = perpuw_max / perpuv_max;
    if (a < 0 || a > 1) || (b < 0 || b > 1)
        intersect_max = 0;
    else
        intersect_max = 1;
        p_max = p + b * v;
    end
else
    intersect_max = 0;
end

if ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
    inv_dir = 1/ray_direction(1);
else
    inv_dir = 1/ray_direction(2);
end

if intersect_max == 1
    if ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
        inv_dir = 1/ray_direction(1);
        t_max = (p_max(1) - ray_origin(1)) * inv_dir;
    else
        inv_dir = 1/ray_direction(2);
        t_max = (p_max(2) - ray_origin(2)) * inv_dir;
    end
end

if intersect_min == 1
    if ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
        inv_dir = 1/ray_direction(1);
        t_min = (p_min(1) - ray_origin(1)) * inv_dir;
    else
        inv_dir = 1/ray_direction(2);
        t_min = (p_min(2) - ray_origin(2)) * inv_dir;
     end
end

if intersect_max == 1 && intersect_min == 0 && ...
        t < t_max && t_max < t_end && ...
        ~approximatelyEqual(t,t_max,1e-12,1e-8)
    tStepTheta = 1;
    tMaxTheta = t_max;
elseif intersect_min == 1 && intersect_max == 0 && ...
        t < t_min && t_min < t_end && ...
        ~approximatelyEqual(t,t_min,1e-12,1e-8)
    tStepTheta = -1;
    tMaxTheta = t_min;
elseif intersect_min == 1 && intersect_max == 1
    % hit the min & max boundary simultaneously.
    if  approximatelyEqual(t_min,t_max,1e-12,1e-8) && ...
            t < t_min && t_min < t_end && ...
            ~approximatelyEqual(t,t_min,1e-12,1e-8)
        tMaxTheta = t_max;
         if verbose
            fprintf("hit origin t\n")
        end
        if ~approximatelyEqual(ray_direction(2),0.0,1e-12,1e-8)
            % change in theta is dependent on the slope of the line
            tStepTheta = - num_angular_sections/2;
        else
            tStepTheta = num_angular_sections/2;
        end
    % hit min bound
    elseif t < t_min && t_min < t_end && ...
            (t_min < t_max || approximatelyEqual(t,t_max,1e-12,1e-8)) && ...
            ~approximatelyEqual(t,t_min,1e-12,1e-8)
        tStepTheta = -1;
        tMaxTheta = t_min;
        if verbose
            fprintf("hit min bound\n")
        end
    % hit max bound
    elseif t < t_max && t_max < t_end && ...
            (t_max < t_min || approximatelyEqual(t,t_min,1e-12,1e-8)) && ...
            ~approximatelyEqual(t,t_max,1e-12,1e-8)
        tStepTheta = 1;
        tMaxTheta = t_max;
        if verbose
            fprintf("hit max bound\n")
         end
     else
        tStepTheta = 0;
        tMaxTheta = t_end;
        if verbose
            fprintf("no hit within time bound\n")
        end
     end
 else
     tStepTheta = 0;
     tMaxTheta = t_end;
     if verbose
        fprintf("no hit within voxel bound\n")
     end
 end


if verbose
    fprintf(['\ntMaxTheta: %d \n' ...
        'tStepTheta: %d \n'], tMaxTheta, tStepTheta);
end
end
