function[current_voxel_ID_angular, p1, P_max_ang] = angular_init(num_angular_sections,...
    sphere_center, sphere_max_radius, current_radius, ray_origin, ray_direction, plane_string)

if strcmp(plane_string, 'xy'); ind1 = 1; ind2 = 2; end
if strcmp(plane_string, 'xz'); ind1 = 1; ind2 = 3; end

[pa] = radial_intersection_points(ray_origin, ray_direction, ...
    sphere_center, current_radius);

delta_ang = 2 * pi/ num_angular_sections;
r = current_radius;
% Create an array of values representing the points of intersection between 
% the lines corresponding to angular voxels boundaries and the initial 
% radial voxel of the ray. Note that spherical coordinates are unnecessary,
% this is just marking the voxel boundaries in the plane. 
i = 1;
k = 0;
trig_ang = zeros(num_angular_sections,2);
while k <= 2*pi 
    trig_ang(i,1) = cos(k);
    trig_ang(i,2) = sin(k);
    k = k + delta_ang;
    i = i + 1;
end
P_max_ang = sphere_max_radius .* trig_ang + [sphere_center(ind1) sphere_center(ind2)];
P_ang = current_radius.* trig_ang + [sphere_center(ind1) sphere_center(ind2)];


% Find the point of intersection between the vector created by the point of
% ray intersection with the initial sphere radius and the sphere center.
ray_origin2d = [ray_origin(ind1), ray_origin(ind2)];
sphere_center2d = [sphere_center(ind1), sphere_center(ind2)];
if approximatelyEqual(ray_origin2d,sphere_center2d,1e-12,1e-8)    
    % If the ray starts at the origin, we need to perturb slightly along its path
    pert_t = 0.1;
    pert_x = ray_origin(ind1) + ray_direction(ind1) * pert_t;
    pert_2 = ray_origin(ind2) + ray_direction(ind2) * pert_t;
    a = sphere_center(ind1) - pert_x;
    b = sphere_center(ind2) - pert_2;
elseif approximatelyEqual(r,sphere_max_radius,1e-12,1e-8)
% If the ray origin is outside the grid, this will snap to the grid and use the
% resulting intersection point as the ray origin for this calculation.
    a = sphere_center(ind1) - pa(ind1);
    b = sphere_center(ind2) - pa(ind2);
else
    a = sphere_center(ind1) - ray_origin(ind1);
    b = sphere_center(ind2) - ray_origin(ind2);
end
l = sqrt(a^2 + b^2);
if approximatelyEqual(l,0.0,1e-12,1e-8)
    % if l is 0.0 then only the z-dir of the ray is non-zero; set the
    % initial point s.t. voxel_ID_theta init is 0.
    p1 = [sphere_center(ind1) + r, sphere_center(ind2)];
else
    p1 = sphere_center2d - (r/l) .* [a b];
end
% This point will lie between two angular voxel boundaries iff the angle between
% it and the angular boundary intersection points along the circle of 
% max radius is obtuse. Equality represents the case when the
% point lies on an angular boundary. 
i = 1;    
while i < length(P_ang)
    d1 = (P_ang(i,1)-p1(1))^2 + (P_ang(i,2)-p1(2))^2;
    d2 = (P_ang(i+1,1)-p1(1))^2 + (P_ang(i+1,2)-p1(2))^2;
    d3 = (P_ang(i,1)-P_ang(i+1,1))^2 + (P_ang(i,2)-P_ang(i+1,2))^2;
    if strictlyLess(d1+d2,d3,1e-12,1e-8) || ...
                approximatelyEqual(d1+d2,d3,1e-12,1e-8)
        current_voxel_ID_angular = i - 1;
        i = length(P_ang);
    end
    i = i + 1;
end

end %end function
