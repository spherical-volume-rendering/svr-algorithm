function[current_voxel_ID_r, r, delta_radius]=radial_init(ray_origin, ray_direction,...
    sphere_center, sphere_max_radius, num_radial_sections, t_begin)

delta_radius = sphere_max_radius / num_radial_sections;

% Determine ray location at t_begin.
p = ray_origin + t_begin .* ray_direction;
p_sphere_vector = sphere_center' - p';
% Find the radial shell containing the ray at t_begin.
r = delta_radius;
while (p_sphere_vector(1)^2 + p_sphere_vector(2)^2 + p_sphere_vector(3)^2 > r^2) && r < sphere_max_radius
    r = r + delta_radius;
end

% If there is a ray/shell intersection, then set the radial voxel ID.
current_voxel_ID_r = 1 + (sphere_max_radius - r) / delta_radius;

end