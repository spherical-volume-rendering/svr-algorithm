function [is_angular_hit, tMaxTheta, tStepTheta] = angular_hit(ray_origin, ray_direction, current_voxel_ID_theta,...
        num_radial_sections, verbose)
% Determines whether an angular hit occurs for the given ray.
% Input:
%    ray_origin: vector of the origin of the ray in cartesian coordinate
%    ray_direction: vector of the direction of the ray in cartesian
%                   coordinate
%    current_voxel_ID: the (angular) ID of current voxel
%    num_radial_sections: number of total radial sections on the grid
% Returns:
%    is_angular_hit: true if an angular crossing has occurred, false otherwise.
%    tMaxTheta: is the time at which a hit occurs for the ray at the next point of intersection.
%    tDeltaTheta: TODO
    if verbose
        fprintf("\nangular_hit.\n");
    end
    
    % First calculate the angular interval that current voxID corresponds
    % to
    delta_theta = 2 * pi / num_radial_sections;
    interval_theta = [current_voxel_ID_theta * delta_theta, (current_voxel_ID_theta + 1) * delta_theta];
    
    % calculate the x and y components that correspond to the angular
    % boundary for the angular interval
    xmin = cos(min(interval_theta));
    xmax = cos(max(interval_theta));
    ymin = sin(min(interval_theta));
    ymax = sin(max(interval_theta));
    
    if verbose
        fprintf("\ninterval_theta: %f \ncurrent_voxel_ID_theta: %f", interval_theta, current_voxel_ID_theta);
        fprintf("\nxmin: %f, xmax: %f \nymin: %f, ymax: %f", xmin, xmax, ymin, ymax);
    end
    
    % solve the systems Az=b to check for intersection
    Amin = [xmin, -ray_direction(1); ymin, -ray_direction(2)];
    Amax = [xmax, -ray_direction(1); ymax, -ray_direction(2)];
    b = transpose([ray_origin(1), ray_origin(2)]);
    zmin = Amin\b; % inv(Amin) * b
    zmax = Amax\b; % inv(Amax) * b
    
    if verbose
        fprintf("\nzmin_r: %f, zmin_t: %f \nzmax_r: %f, zmax_t: %f", zmin(1), zmin(2), zmax(1), zmax(2));    
    end
    
    
    % We need the radius (r = z[0]) and time (t = z[0]) to be positive or
    % else the intersection is null
    is_angular_hit = true;
    if zmin(1) < 0 || zmin(2) < 0
        is_angular_hit = false;
        tMaxTheta = -inf;
        tStepTheta = -inf;
        return;
    end
    if zmax(1) < 0 || zmin(2) < 0
        is_angular_hit = false;
                tMaxTheta = -inf;
        tStepTheta = -inf;
        return;
    end
    
    % If we hit the min boundary then we decrement theta, else increment;
    % assign tmaxtheta
    if zmin(1) < 0 || zmin(2) < 0 
        tStepTheta = -1;
        tMaxTheta = zmin(1);
    else
        tStepTheta = 1;
        tMaxTheta = zmax(1);
    end
    if verbose
        fprintf([...
                 'tMaxTheta: %d \n' ...
                 'is_angular_hit: %d \n' ...
                 'tStepTheta: %d \n'], tMaxTheta, is_angular_hit, tStepTheta);
    end
end
