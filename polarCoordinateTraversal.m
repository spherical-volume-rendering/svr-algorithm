function polarCoordinateTraversal(min_bound, max_bound, ray_origin, ray_direction, circle_center, circle_max_radius, num_radial_sections, num_angular_sections, t_begin, t_end, verbose)
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
            quiver(ray_origin_x, ray_origin_y, ray_direction_x, ray_direction_y, t_begin - 0.0);
        end
        
        text(ray_start_x, ray_start_y, ' ray start');
        text(ray_end_x, ray_end_y, ' ray end');
        plot(ray_end_x, ray_end_y, 'k.', 'MarkerSize', 10);
        plot(ray_start_x, ray_start_y, 'k.', 'MarkerSize', 10);
        quiver(ray_start_x, ray_start_y, ray_direction_x, ray_direction_y, t_end - t_begin);
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
        % TODO(cpg49)
    end
    
    % If the ray does intersect the grid, calculate the 
    % voxel_ID(r, theta) of entry. Note that the outermost shell is the
    % first radial section. We are then guaranteed:
    % 1 <= voxel_ID_r <= num_radial_sections
    voxel_ID_r = 1;
    
    % Similarly, we can define unique angular sections. This guarantees us:
    % 1 <= voxel_ID_theta <= num_angular_sections
    voxel_ID_theta = floor(atan2(ray_start_y, ray_start_x) * num_angular_sections / (2 * pi));
    if voxel_ID_theta < 0
        voxel_ID_theta = num_angular_sections + voxel_ID_theta;
    end
    
    t = t_hit;
    
    
    while t < t_exit
        % [is_angular_hit, tMax_theta, delta_theta] = angular_hit(ray_origin, ray_direction, voxel_ID_theta);
        break;
        % TODO
    end
    
end

function [is_radial_hit, t_hit, t_exit] = radial_hit(ray_start, ray_direction, max_radius)
% Returns is_radial_hit = true if a radial_hit has occurred, false otherwise.
% t_hit is the location where the ray has entered.
% t_exit is the location where the ray exits.
    assert(false);
    % TODO
end

function [is_angular_hit, tMax_theta, delta_theta] = angular_hit(ray_start, ray_direction, voxel_ID_theta)
    % Determines whether the ray hits an angular section. Note that
    % delta_theta = +1 or -1; either we increase or decrease in angular
    % direction.
    assert(false);
    % TODO
    
end
