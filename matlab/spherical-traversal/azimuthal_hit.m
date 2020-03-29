function [tMaxPhi, tStepPhi] = azimuthal_hit(ray_origin, ray_direction, ...
        current_voxel_ID_phi, num_azimuthal_sections, sphere_center, t, ...
        verbose)
% Determines whether an azimuthal hit occurs for the given ray.
% The azimuthal voxels lie on the X-Z plane.
% Input:
%    ray_origin: vector of the origin of the ray in cartesian coordinates.
%    ray_direction: vector of the direction of the ray in cartesian coordinates.
%    current_voxel_ID_phi: the (azimuthal) ID of current voxel.
%    num_azimuthal_sections: number of total azimuthal sections of the grid.
%    sphere_center: The center of the sphere.
%    t: The current time parameter of the ray.
%    verbose: Determines whether debug print statements are enabled.
%
% Returns:
%    tMaxPhi: is the time at which a hit occurs for the ray at the next point of intersection.
%    tStepPhi: The direction the azimuthal voxel steps. +1, -1, or 0.
if verbose
    fprintf("\n-- angular_hit --")
end


% First calculate the angular interval that current voxelID corresponds to.
delta_phi = 2 * pi / num_azimuthal_sections;
interval_phi = [current_voxel_ID_phi * delta_phi, (current_voxel_ID_phi + 1) * delta_phi];

if verbose
    fprintf("\nCurrent Voxel ID Theta: %d\n", current_voxel_ID_phi)
end


% Calculate the x and y components that correspond to the angular
% boundary for the angular interval.
xmin = cos(min(interval_phi));
xmax = cos(max(interval_phi));
zmin = sin(min(interval_phi));
zmax = sin(max(interval_phi));
tol = 10^-12;
% Solve the systems Au=b to check for intersection.
Amin = [xmin, -ray_direction(1); zmin, -ray_direction(3)];
Amax = [xmax, -ray_direction(1); zmax, -ray_direction(3)];
b = [ray_origin(1)-sphere_center(1), ray_origin(3)-sphere_center(3)]';
if abs(det(Amin)) < tol
    umax = Amax \ b;
    umin = [0 ; 0];
elseif abs(det(Amax)) < tol
    umin = Amin\b;
    umax = [0 ; 0];
else
    umin = Amin\b;
    umax = Amax\b;
end

if umin(1) > tol && umin(2) > t
    % If we hit the min boundary then we decrement theta, else increment;
    % assign tMaxPhi.
    tStepPhi = -1;
    tMaxPhi = umin(2);
    if verbose
        fprintf("hit min bound\n")
    end
elseif umax(1) > tol && umax(2) > t
    tStepPhi = 1;
    tMaxPhi = umax(2);
    if verbose
        fprintf("hit max bound\n")
    end
elseif abs(umax(1) - umin(1)) < tol && umax(2) -t > tol
        % hit the min & max boundary simultaneously. 
        tMaxPhi = umax(2);
        if verbose
          fprintf("hit origin max t\n")
        end
    if ray_direction(2) > tol
        % change in theta is dependent on the slope of the line
        interval_phi = interval_phi - pi;
        tStepPhi = - abs(current_voxel_ID_phi - interval_phi(1)/delta_phi);
    else
        interval_phi = interval_phi + pi;
        tStepPhi = abs(current_voxel_ID_phi - interval_phi(1)/delta_phi);
    end
elseif abs(umax(1) - umin(1)) < tol && umin(2) -t > tol
    % hit the min & max boundary simultaneously. 
    tMaxPhi = umin(2);
    if verbose
        fprintf("hit origin min t\n")
    end
    if ray_direction(2) > tol
        % change in theta is dependent on the slope of the line
        interval_phi = interval_phi - pi;
        tStepPhi = - abs(current_voxel_ID_phi - interval_phi(1)/delta_phi);
    else
        interval_phi = interval_phi + pi;
        tStepPhi = abs(current_voxel_ID_phi - interval_phi(1)/delta_phi);
    end
else
    tStepPhi = 0;
    tMaxPhi = inf;
    if verbose
        fprintf("no hit\n")
    end
end
end
