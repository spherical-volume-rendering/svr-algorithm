function[p1,p2,discr]=radial_intersection_points(ray_origin, ray_direction,...
        sphere_center, radius)
    
    ray_sphere_vector = sphere_center' - ray_origin' ;
    r = radius; 
    ray_unit_vector = 1 / sqrt(ray_direction(1)^2 + ray_direction(2)^2 + ray_direction(3)^2)...
       .* [ray_direction(1);  ray_direction(2); ray_direction(3)]';
    v = dot(ray_sphere_vector, ray_unit_vector);
    discr = r^2 - (dot(ray_sphere_vector,ray_sphere_vector) - v^2);
    d = sqrt(discr);
    p1 = ray_origin + (v-d) .* ray_unit_vector;
    p2 = ray_origin + (v+d) .* ray_unit_vector;
    
end