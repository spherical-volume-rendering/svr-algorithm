function polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, circle_center, ...
        circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose)
% Input:
%    min_bound: The lower left corner of the bounding box.
%    max_bound: The upper right corner of the bounding box.
%    ray origin: The origin of the ray in (x, y) coordinates.
%    ray direction: The direction of the ray in (x, y) coordinates.
%    circle_center: The x, y location of the center of the circle.
%    circle_max_radius: The largest that encompasses the circle.
%    num_radial_sections: The number of radial sections in the circle.
%    num_angular_sections: The number of angular sections in the circle.
%    t_begin: The beginning time of the ray.
%    t_end: The end time of the ray.
%
% Requires:
%    max_bound > min_bound
%    circle center is within max_bound and min_bound.
%    t_end > t_begin >= 0.0
%    circle_max_radius > 0
%    num_radial_sections > 0
%    num_angular_sections > 0
%
% Notes: 
%    Currently under construction.
    close all;
    circle_center_x = circle_center(1);
    circle_center_y = circle_center(2);
    ray_origin_x = ray_origin(1);
    ray_origin_y = ray_origin(2);
    ray_direction_x = ray_direction(1);
    ray_direction_y = ray_direction(2);
    
    min_bound_x = min_bound(1);
    min_bound_y = min_bound(2);
    max_bound_x = max_bound(1);
    max_bound_y = max_bound(2);
    
    ray_start = ray_origin + t_begin * ray_direction;
    ray_start_x = ray_start(1);
    ray_start_y = ray_start(2);
    
    ray_end = ray_origin + t_end * ray_direction;
    ray_end_x = ray_end(1);
    ray_end_y = ray_end(2);
    
    if (verbose)
        figure;
        hold on;
        title('Polar Coordinate Voxel Traversal')
        
        if (t_begin ~= 0.0)
            % Mark the ray origin if the time does not start at 0.0
            text(ray_origin_x, ray_origin_y, ' ray origin');
            plot(ray_origin_x, ray_origin_y, 'k.', 'MarkerSize', 10);
            quiver(ray_origin_x, ray_origin_y, ray_direction_x, ray_direction_y, t_begin - 0.0, 'LineWidth', 1.5);
        end
        
        % Draw the ray.
        text(ray_start_x, ray_start_y, ' ray start');
        text(ray_end_x, ray_end_y, ' ray end');
        plot(ray_end_x, ray_end_y, 'k.', 'MarkerSize', 10);
        plot(ray_start_x, ray_start_y, 'k.', 'MarkerSize', 10);
        quiver(ray_start_x, ray_start_y, ray_direction_x, ray_direction_y, t_end - t_begin, 'LineWidth', 1.5);
        
        % Draw the axis.
        axis tight;
        xlim([min_bound_x, max_bound_x]);
        ylim([min_bound_y, max_bound_y]);
        xlabel('x');
        ylabel('y');
        grid on;
        
        % Draw the radial sections.
        current_max_radius = circle_max_radius;
        delta_radius = circle_max_radius / num_radial_sections;
        for k = 1:num_radial_sections
            viscircles(circle_center, current_max_radius, 'LineStyle', '--', 'Color', '#7E2F8E', 'LineWidth', 1);
            current_max_radius = current_max_radius - delta_radius;
        end
        
        % Draw the angular sections.
        N = num_angular_sections;
        section = 2 * pi / num_angular_sections;
        for ii = 1:N
              t = linspace(section * (ii - 1), section * (ii));
              x = circle_max_radius*cos(t) + circle_center_x;
              y = circle_max_radius*sin(t) + circle_center_y;
              x = [x circle_center_x x(1)];
              y = [y circle_center_y y(1)];
              line(x, y, 'LineStyle', '--', 'Color', '#7E2F8E', 'LineWidth', 0.5);
        end
    end
    
    % INITIALIZATION PHASE
    %  I. Calculate Voxel ID R.
    delta_radius = circle_max_radius / num_radial_sections;
    current_position = (ray_start_x - circle_center_x)^2 + (ray_start_y - circle_center_y)^2;
    if current_position > circle_max_radius
        voxel_ID_r = 1;
    else
        current_delta_radius = delta_radius;
        current_voxel_ID_r = num_radial_sections;
        while (current_position < current_delta_radius)
            current_voxel_ID_r = current_voxel_ID_r - 1;
            current_delta_radius = current_delta_radius + delta_radius;
        end
    end
    
    % II. Calculate Voxel ID Theta.
    current_voxel_ID_theta = floor(atan2(ray_start_y, ray_start_x) * num_angular_sections / (2 * pi));
    if current_voxel_ID_theta < 0
        current_voxel_ID_theta = num_angular_sections + current_voxel_ID_theta;
    end
    
    % TRAVERSAL PHASE
    current_ray_position = ray_start;
    t = t_begin;
    while t < t_end
        % Calculate tMaxR (using radial_hit) 
        [is_radial_hit, tMaxR, tDeltaR, new_voxel_ID_r] = radial_hit(current_ray_position, ray_direction, ...
            current_voxel_ID_r, circle_center, circle_max_radius, delta_radius);
        
        % TODO: Calculate tMaxTheta (using angular_hit)
        
        % Compare tMaxTheta, tMaxR
        if tMaxTheta < tMaxR
            t = t + tDeltaTheta;
        else
            t = t + tDeltaR;
        end
        % Update new Voxel IDs.
        current_voxel_ID_r = new_voxel_ID_r;
        % TODO: update current theta voxel.
        
    end
    
end

function [is_radial_hit, tMaxR, tDeltaR, new_voxel_ID_r] = radial_hit(current_ray_position, ray_direction, ...
        current_radial_voxel, circle_center, circle_max_radius, delta_radius)
% Determines whether a radial hit occurs for the given ray.
% Input:
%    current_ray_position: The first location of the ray within the circle.
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
%    tDeltaR: The different between the ray new position and the ray's
%    initial position along its direction.
%    new_voxel_ID_r: The new voxel ID that the ray is located in. If the
%    ray hasn't changed, this remains the old voxel ID.
    ray_pos_x = current_ray_position(1);
    ray_pos_y = current_ray_position(2);
    ray_direction_x = ray_direction(1);
    ray_direction_y = ray_direction(2);
    circle_center_x = circle_center(1);
    circle_center_y = circle_center(2);
    current_radius = circle_max_radius - (delta_radius * (current_radial_voxel - 1));
    
    % (1)   (x - circle_center_x)^2 + (y - circle_center_y)^2 = current_radius^2
    % (2)    x = ray_origin_x + ray_direction_x(t)
    % (3)    y = ray_origin_y + ray_direction_y(t)
    % Plug in x, y in equation (1), then solve for t.
    % To get point of intersection, plug t back in parametric equation of a ray.
    tMaxR = solve((ray_pos_x + ray_direction_x * t - circle_center_x)^2 + ...
        (ray_pos_y + ray_direction_y * t - circle_center_y)^2 ...
        - current_radius^2 == 0, t);
    new_x_position = ray_pos_x + ray_direction_x * tMaxR;
    new_y_position = ray_pos_y + ray_direction_y * tMaxR;
    
    % Determine whether is has switched to a new radial voxel.
    current_position = (new_x_position - circle_center_x)^2 + (new_y_position - circle_center_y)^2;
    if current_position >= current_radius + delta_radius
        new_voxel_ID_r = current_radial_voxel + 1;
        is_radial_hit = true;
    elseif current_position < current_radius
            new_voxel_ID_r = current_radial_voxel - 1;
            is_radial_hit = true;
    else
        new_voxel_ID_r = current_radial_voxel;
        is_radial_hit = false;
    end
    
    % Calculate tDeltaR.
    tDeltaR = 0;
    if is_radial_hit
        tDeltaR = current_ray_position - [new_x_position new_y_position];
    end
end

function [is_angular_hit, tMaxTheta, tDeltaTheta] = angular_hit(current_ray_pos, ray_direction, current_voxel_ID_theta)
% Determines whether an angular hit occurs for the given ray.
% UPDATE DOCUMENTATION.
% Returns:
%    is_angular_hit: true if an angular crossing has occurred, false otherwise.
%    tMaxTheta: is the time at which a hit occurs for the ray at the next point of intersection.
%    tDeltaTheta: TODO

    % TODO: Implement
    assert(false); 
    % 1. Parametric equation of ray using origin and direction.
    % Using voxID, calculate the theta interval. 
    % i.e. voxID_theta = voxIDtheta * num_angular_sections / 2pi or
    % something like this
    %
    % Calculate x, y using rcos(max[voxID_theta]), rsin(max[voxID_theta])
    
end
