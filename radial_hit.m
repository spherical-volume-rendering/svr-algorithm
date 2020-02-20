function [is_radial_hit, tMaxR, tStepR] = ...
        radial_hit(ray_origin, ray_direction, ...
        current_radial_voxel, circle_center, ...
        circle_max_radius, delta_radius, verbose)
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
%    tStepR: The direction of step into the next radial voxel, 0, +1, -1
%    new_ray_position: The (x,y) coordinate of the ray after the traversal.
    ray_direction_x = ray_direction(1);
    ray_direction_y = ray_direction(2);
    circle_center_x = circle_center(1);
    circle_center_y = circle_center(2);
    ray_origin_x = ray_origin(1);
    ray_origin_y = ray_origin(2);
    current_radius = circle_max_radius - (delta_radius * (current_radial_voxel - 1));
    
    if verbose
        fprintf('\nradial_hit. \nCurrent Radial Voxel: %d\n', current_radial_voxel);
    end
        
    % (1)   (x - circle_center_x)^2 + (y - circle_center_y)^2 = current_radius^2
    % (2)    x = ray_origin_x + ray_direction_x(t)
    % (3)    y = ray_origin_y + ray_direction_y(t)
    % Plug in x, y in equation (1), then solve for t.
    % To get point of intersection, plug t back in parametric equation of a ray.
    syms cT; % current time

    intersections_t = solve((ray_origin_x + ray_direction_x * cT - circle_center_x)^2 + ...
        (ray_origin_y + ray_direction_y * cT - circle_center_y)^2 ...
        - (current_radius - delta_radius)^2 == 0, cT);
    
    if isempty(double(subs(intersections_t)))
        % Case where no intersection between ray and new radial voxel
        % occurs
        is_radial_hit = false;
        tStepR = 0;
    end
    
    tMaxR = min(double(subs(intersections_t)));
    
    new_x_position = ray_origin_x + ray_direction_x * tMaxR;
    new_y_position = ray_origin_y + ray_direction_y * tMaxR;
    
    if verbose
        fprintf('tMaxR %f\n', tMaxR);
        fprintf('New position: (%f, %f). Calculated by ray_origin + ray_direction * tMaxR.\n ', new_x_position, new_y_position);
    end
    
    distance_from_origin = (new_x_position - circle_center_x)^2 + (new_y_position - circle_center_y);
    
    if verbose
        fprintf("Distance from origin ^2: %f\n",  distance_from_origin);
        fprintf("Current radius ^2: %f\n", current_radius^2);
        fprintf("(Current radius - delta_radius) ^2: %f\n", (current_radius - delta_radius)^2);
    end
    
    if  distance_from_origin >= current_radius^2
        is_radial_hit = true;
        tStepR = -1;
        if verbose
            text(new_x_position, new_y_position, 'POI_r');
            fprintf('Ray moving toward voxel closer to perimeter (outward).\n');
        end
    else  
        % distance_from_origin < (current_radius - delta_radius)^2
        is_radial_hit = true;
        tStepR = +1;
        if verbose
            text(new_x_position, new_y_position, 'POI_r');
            fprintf('Ray moving toward voxel closer to center (inward).\n');
        end
    end
    
    if verbose
        fprintf(['new_voxel_ID_r: %d \n' ...
            'is_radial_hit: %d \n' ...
            'tStepR: %d \n'], current_radial_voxel + tStepR, is_radial_hit, tStepR);
    end
end
