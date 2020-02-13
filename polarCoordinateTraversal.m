function polarCoordinateTraversal(ray_origin, ray_direction, max_radius, num_radial_sections, num_angular_sections)
% Input:
%    ray origin: The origin of the ray in (x, y) coordinates.
%    ray direction: The direction of the ray in (x, y) coordinates.
%    max_radius: The largest that encompasses the circle.
%    num_radial_sections: The number of radial sections in the circle.
%    num_angular_sections: The number of angular sections in the circle.
%
% Requires:
%    max_radius > 0
%    num_radial_sections > 0
%    num_angular_sections > 0
%
% Notes: 
%    Currently under construction.

    ray_origin_x = ray_origin(1);
    ray_origin_y = ray_origin(2);
    ray_direction_x = ray_direction(1);
    ray_direction_y = ray_direction(2);
    
    [is_radial_hit, t_hit, t_exit] = radial_hit(ray_origin, ray_direction, max_radius);
    
    if ~is_radial_hit
        % The ray does not hit the first radius.
        return;
    end
    
    % If the ray does intersect the grid, calculate the 
    % voxel_ID(r, theta) of entry. Note that the outermost shell is the
    % first radial section. We are then guaranteed:
    % 1 <= voxel_ID_r <= num_radial_sections
    voxel_ID_r = 1;
    
    % Similarly, we can define unique angular sections. This guarantees us:
    % 1 <= voxel_ID_theta <= num_angular_sections
    voxel_ID_theta = floor(atan2(ray_origin_y, ray_origin_x) * num_angular_sections / (2 * pi));
    if voxel_ID_theta < 0
        voxel_ID_theta = num_angular_sections + voxel_ID_theta;
    end
    
    t = t_hit;
    
    
    while t < t_exit
        [is_angular_hit, tMax_theta, delta_theta] = angular_hit(ray_origin, ray_direction, voxel_ID_theta);
        
        % TODO
    end
    
end

function [is_radial_hit, t_hit, t_exit] = radial_hit(ray_origin, ray_direction, max_radius)
% Returns is_radial_hit = true if a radial_hit has occurred, false otherwise.
% t_hit is the location where the ray has entered.
% t_exit is the location where the ray exits.
    assert(false);
    % TODO
end

function [is_angular_hit, tMax_theta, delta_theta] = angular_hit(ray_origin, ray_direction, voxel_ID_theta)
    % Determines whether the ray hits an angular section. Note that
    % delta_theta = +1 or -1; either we increase or decrease in angular
    % direction.
    assert(false);
    % TODO
    
end
