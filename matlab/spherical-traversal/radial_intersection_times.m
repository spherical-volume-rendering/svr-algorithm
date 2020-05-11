function[t1,t2]=radial_intersection_times(ray_origin, ray_direction,...
        sphere_center, radius)
    
    [pa, pb] = radial_intersection_points(ray_origin, ray_direction,...
        sphere_center, radius);
    % Calculate the time of entrance and exit of the ray.
    if ~approximatelyEqual(ray_direction(2),0.0,1e-12,1e-8)
        % Use the y-direction if it is non-zero.
        t1 = (pa(2) - ray_origin(2)) / ray_direction(2);
        t2 = (pb(2) - ray_origin(2)) / ray_direction(2);
    elseif ~approximatelyEqual(ray_direction(1),0.0,1e-12,1e-8)
        % Use the x-direction if it is non-zero.
        t1 = (pa(1) - ray_origin(1)) / ray_direction(1);
        t2 = (pb(1) - ray_origin(1)) / ray_direction(1);
    else
        % Use the z-direction if it is non-zero.
        t1 = (pa(3) - ray_origin(3)) / ray_direction(3);
        t2 = (pb(3) - ray_origin(3)) / ray_direction(3);
    end
    
end