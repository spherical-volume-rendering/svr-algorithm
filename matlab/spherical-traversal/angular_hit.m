function [tMaxTheta, tStepTheta] = angular_hit(ray_origin, ray_direction, current_voxel_ID_theta,...
    sphere_center, P_max, sphere_max_radius, t, t_end, verbose)
% Determines whether an angular hit occurs for the given ray.
% Input:
%    ray_origin: vector of the origin of the ray in cartesian coordinate.
%    ray_direction: vector of the direction of the ray in cartesian coordinate.
%    current_voxel_ID_theta: the (angular) ID of current voxel.
%    num_angular_sections: number of total angular sections on the grid.
%    sphere_center: The center of the circle.
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
pointArray = [P_max(current_voxel_ID_theta+1,:) ; P_max(current_voxel_ID_theta+2,:)];
% Get current ray location to use as start point
p = ray_origin + ray_direction * t;
% Get ray exit point to use as end point
p_f = ray_origin + ray_direction * t_end;
% Caculate the ray segment vector
v = p_f - p;
% Calculate the voxel boundary vectors
u_min = [sphere_center(1) sphere_center(2)] - pointArray(1,:);
u_max = [sphere_center(1) sphere_center(2)] - pointArray(2,:);
% http://geomalgorithms.com/a05-_intersect-1.html#intersect2D_2Segments()
w_min = pointArray(1,:) - [p(1) p(2)];
w_max = pointArray(2,:) - [p(1) p(2)];
perpuv_min = u_min(1) * v(2) - u_min(2) * v(1);
perpuv_max = u_max(1) * v(2) - u_max(2) * v(1);
perpuw_min = u_min(1) * w_min(2) - u_min(2) * w_min(1);
perpuw_max = u_max(1) * w_max(2) - u_max(2) * w_max(1);
perpvw_min = v(1) * w_min(2) - v(2) * w_min(1);
perpvw_max = v(1) * w_max(2) - v(2) * w_max(1);

% check to see if ray is parallel to each angular voxel boundary
if approximatelyEqual(perpuv_min,0.0,1e-12,1e-8) 
    parallel_min = 1;
    % check to see if the ray is colliner with the angular voxel boundary
    if approximatelyEqual(perpuw_min,0.0,1e-12,1e-8) && ...
            approximatelyEqual(perpvw_min,0.0,1e-12,1e-8)
        collinear_min = 1;
        p_min = sphere_center;
    else
        collinear_min = 0;
    end
else
    parallel_min = 0;
    collinear_min = 0;
end
if approximatelyEqual(perpuv_max,0.0,1e-12,1e-8)
    parallel_max = 1;
    if approximatelyEqual(perpuw_max,0.0,1e-12,1e-8) && ...
            approximatelyEqual(perpvw_max,0.0,1e-12,1e-8)
        collinear_max = 1;
        p_max = sphere_center;
    else
        collinear_max = 0;
    end
else
    parallel_max = 0;
    collinear_max = 0;
end
% Determine intersection points, p_min and p_max
if parallel_min == 0
    a = perpvw_min / perpuv_min;
    b = perpuw_min / perpuv_min;
    if (strictlyLess(a,0.0,1e-12,1e-8) || strictlyLess(1.0,a,1e-12,1e-8)) || ...
            (strictlyLess(b,0.0,1e-12,1e-8) || strictlyLess(1.0,b,1e-12,1e-8))
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
    if (strictlyLess(a,0.0,1e-12,1e-8) || strictlyLess(1.0,a,1e-12,1e-8)) || ...
            (strictlyLess(b,0.0,1e-12,1e-8) || strictlyLess(1.0,b,1e-12,1e-8))
        intersect_max = 0;
    else
        intersect_max = 1;
        p_max = p + b * v;
    end
else
    intersect_max = 0;
end
% Compute inverses of slopes
if ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
    inv_dir = 1/ray_direction(1);
else
    inv_dir = 1/ray_direction(2);
end
% Compute intersection times
if intersect_max == 1 || collinear_max == 1
    if ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
        t_max = (p_max(1) - ray_origin(1)) * inv_dir;
    else
        t_max = (p_max(2) - ray_origin(2)) * inv_dir;
    end
end
if intersect_min == 1 || collinear_min == 1
    if ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
        t_min = (p_min(1) - ray_origin(1)) * inv_dir;
    else
        t_min = (p_min(2) - ray_origin(2)) * inv_dir;
     end
end
% Determine tStepTheta, tMaxTheta
if intersect_max == 1 && intersect_min == 0 && collinear_min == 0 &&...
        strictlyLess(t,t_max,1e-12,1e-8) && strictlyLess(t_max,t_end,1e-12,1e-8)
    tStepTheta = 1;
    tMaxTheta = t_max;
elseif intersect_min == 1 && intersect_max == 0 && collinear_max == 0 && ...
           strictlyLess(t,t_min,1e-12,1e-8) && ...
           strictlyLess(t_min,t_end,1e-12,1e-8)
    tStepTheta = -1;
    tMaxTheta = t_min;
elseif (intersect_min == 1 && intersect_max == 1) || ...
        (intersect_min == 1 && collinear_max == 1) || ...
        (collinear_min == 1 && intersect_max == 1)
    % hit the min & max boundary simultaneously.
    if  approximatelyEqual(t_min,t_max,1e-12,1e-8) && ...
           strictlyLess(t,t_min,1e-12,1e-8) && ...
            strictlyLess(t_min,t_end,1e-12,1e-8)
        tMaxTheta = t_max;
         if verbose
            fprintf("hit origin t\n")
         end
         pert_t = 0.1;
         pert_x = sphere_center(1) + ray_direction(1) * pert_t;
         pert_y = sphere_center(2) + ray_direction(2) * pert_t;
         a = sphere_center(1) - pert_x;
         b = sphere_center(2) - pert_y;
         l = sqrt(a^2 + b^2);
         sphere_center_xy = [sphere_center(1), sphere_center(2)];
         p1 = sphere_center_xy - (sphere_max_radius/l) .* [a b];
         i = 1;    
        while i < length(P_max)
            d1 = (P_max(i,1)-p1(1))^2 + (P_max(i,2)-p1(2))^2;
            d2 = (P_max(i+1,1)-p1(1))^2 + (P_max(i+1,2)-p1(2))^2;
            d3 = (P_max(i,1)-P_max(i+1,1))^2 + (P_max(i,2)-P_max(i+1,2))^2;
            if strictlyLess(d1+d2,d3,1e-12,1e-8) || ...
                approximatelyEqual(d1+d2,d3,1e-12,1e-8)
                next_vox = i - 1 ;
                i = length(P_max);
            end
            i = i + 1;
        end
        step = abs(current_voxel_ID_theta - next_vox);
        if strictlyLess(ray_direction(2),0.0,1e-12,1e-8)
            tStepTheta = step;
        elseif strictlyLess(0.0,ray_direction(2),1e-12,1e-8)
            tStepTheta = -step;
        else
            if strictlyLess(ray_direction(1),0.0,1e-12,1e-8)
                tStepTheta = step;
            else
                tStepTheta = -step;
            end
        end
    % hit min bound
    elseif strictlyLess(t,t_min,1e-12,1e-8) && ...
           strictlyLess(t_min,t_end,1e-12,1e-8) && ...
            (strictlyLess(t_min,t_max,1e-12,1e-8) || ...
            approximatelyEqual(t,t_max,1e-12,1e-8))
        tStepTheta = -1;
        tMaxTheta = t_min;
        if verbose
            fprintf("hit min bound\n")
        end
    % hit max bound
   elseif strictlyLess(t,t_max,1e-12,1e-8) && ...
           strictlyLess(t_max,t_end,1e-12,1e-8) && ...
            (strictlyLess(t_max,t_min,1e-12,1e-8) || ...
            approximatelyEqual(t,t_min,1e-12,1e-8))
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
