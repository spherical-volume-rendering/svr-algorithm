function [tMaxPhi, tStepPhi] = azimuthal_hit(ray_origin, ray_direction, ...
        current_voxel_ID_phi, num_angular_sections, sphere_center, t, ...
        verbose)
% Determines whether an azimuthal hit occurs for the given ray.
% The azimuthal voxels lie on the X-Z plane.
% Input:
%    ray_origin: vector of the origin of the ray in cartesian coordinates.
%    ray_direction: vector of the direction of the ray in cartesian coordinates.
%    current_voxel_ID_phi: the (azimuthal) ID of current voxel.
%    num_angular_sections: number of total angular sections on the grid.
%    sphere_center: The center of the sphere.
%    t: The current time parameter of the ray.
%    verbose: Determines whether debug print statements are enabled.
%
% Returns:
%    tMaxPhi: is the time at which a hit occurs for the ray at the next point of intersection.
%    tStepPhi: The direction the azimuthal voxel steps. +1, -1, or 0.
end
