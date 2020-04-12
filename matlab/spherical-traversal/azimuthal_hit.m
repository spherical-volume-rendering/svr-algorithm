function [tMaxPhi, tStepPhi] = azimuthal_hit(ray_origin, ray_direction, current_voxel_ID_phi,...
    num_azimuthal_sections, sphere_center, pointArray, t, t_end, verbose)
% Determines whether an azimuthal hit occurs for the given ray.
% Input:
%    ray_origin: vector of the origin of the ray in cartesian coordinate.
%    ray_direction: vector of the direction of the ray in cartesian coordinate.
%    current_voxel_ID_phi: the (azimuthal) ID of current voxel.
%    num_azimuthal_sections: number of total azimuthal sections on the grid.
%    sphere_center: The center of the circle.
%    pointArray:
%    t: The current time parameter of the ray.
%    t_end:
%    verbose: Determines whether debug print statements are enabled.
%
% Returns:
%    tMaxPhi: is the time at which a hit occurs for the ray at the next point of intersection.
%    tStepPhi: The direction the phi voxel steps. +1, -1, or 0.
if verbose
    fprintf("\n-- azimuthal_hit --")
end

if verbose
    fprintf("\nCurrent Voxel ID Phi: %d\n", current_voxel_ID_phi)
end

% Get current ray location to use as start point
p = ray_origin + ray_direction * t;
% Get ray exit point to use as end point
p_f = ray_origin + ray_direction * t_end;
% Caculate the ray segment vector
v = p_f - p;
% Calculate the voxel boundary vectors
u_min = [sphere_center(1) sphere_center(3)] - pointArray(1,:);
u_max = [sphere_center(1) sphere_center(3)] - pointArray(2,:);
% http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
w_min = pointArray(1,:) - [p(1) p(3)];
w_max = pointArray(2,:) - [p(1) p(3)];
perpuv_min = u_min(1) * v(2) - u_min(2) * v(1);
perpuv_max = u_max(1) * v(2) - u_max(2) * v(1);
perpuw_min = u_min(1) * w_min(2) - u_min(2) * w_min(1);
perpuw_max = u_max(1) * w_max(2) - u_max(2) * w_max(1);
perpvw_min = v(1) * w_min(2) - v(2) * w_min(1);
perpvw_max = v(1) * w_max(2) - v(2) * w_max(1);

% check to see if ray is parallel to each azimuthal voxel boundary
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
    inv_dir = 1/ray_direction(3);
end

if intersect_max == 1
    if ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
        inv_dir = 1/ray_direction(1);
        t_max = (p_max(1) - ray_origin(1)) * inv_dir;
    else
        inv_dir = 1/ray_direction(3);
        t_max = (p_max(2) - ray_origin(3)) * inv_dir;
    end
end

if intersect_min == 1
    if ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
        inv_dir = 1/ray_direction(1);
        t_min = (p_min(1) - ray_origin(1)) * inv_dir;
    else
        inv_dir = 1/ray_direction(3);
        t_min = (p_min(2) - ray_origin(3)) * inv_dir;
     end
end

if intersect_max == 1 && intersect_min == 0 && ...
        t < t_max && t_max < t_end && ...
        ~approximatelyEqual(t,t_max,1e-12,1e-8)
    tStepPhi = 1;
    tMaxPhi = t_max;
elseif intersect_min == 1 && intersect_max == 0 && ...
        t < t_min && t_min < t_end && ...
        ~approximatelyEqual(t,t_min,1e-12,1e-8)
    tStepPhi = -1;
    tMaxPhi = t_min;
elseif intersect_min == 1 && intersect_max == 1
    % hit the min & max boundary simultaneously.
    if  approximatelyEqual(t_min,t_max,1e-12,1e-8) && ...
            t < t_min && t_min < t_end && ...
            ~approximatelyEqual(t,t_min,1e-12,1e-8)
        tMaxPhi = t_max;
         if verbose
            fprintf("hit origin t\n")
        end
        if ~approximatelyEqual(ray_direction(3),0.0,1e-12,1e-8)
            % change in phi is dependent on the slope of the line
            tStepPhi = - num_azimuthal_sections/2;
        else
            tStepPhi = num_azimuthal_sections/2;
        end
    % hit min bound
    elseif t < t_min && t_min < t_end && ...
            (t_min < t_max || approximatelyEqual(t,t_max,1e-12,1e-8)) && ...
            ~approximatelyEqual(t,t_min,1e-12,1e-8)
        tStepPhi = -1;
        tMaxPhi = t_min;
        if verbose
            fprintf("hit min bound\n")
        end
    % hit max bound
    elseif t < t_max && t_max < t_end && ...
            (t_max < t_min || approximatelyEqual(t,t_min,1e-12,1e-8)) && ...
            ~approximatelyEqual(t,t_max,1e-12,1e-8)
        tStepPhi = 1;
        tMaxPhi = t_max;
        if verbose
            fprintf("hit max bound\n")
         end
     else
        tStepPhi = 0;
        tMaxPhi = t_end;
        if verbose
            fprintf("no hit within time bound\n")
        end
     end
 else
     tStepPhi = 0;
     tMaxPhi = t_end;
     if verbose
        fprintf("no hit within voxel bound\n")
     end
 end


if verbose
    fprintf(['\ntMaxPhi: %d \n' ...
        'tStepPhi: %d \n'], tMaxPhi, tStepPhi);
end
end
