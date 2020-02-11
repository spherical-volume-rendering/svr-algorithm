function [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax)
%  Ray/box intersection using the Smits' algorithm
%
% Input:
%    origin.
%    direction.
%    box = (vmin,vmax)
% Output:
%    flag: (0) Reject, (1) Intersect.
%    tmin: distance from the ray origin.
% Author: 
%    Jesus Mena

ray_direction_x = direction(1);
ray_direction_y = direction(2);
ray_direction_z = direction(3);

ray_origin_x = origin(1);
ray_origin_y = origin(2);
ray_origin_z = origin(3);

minBound_x = vmin(1);
minBound_y = vmin(2);
minBound_z = vmin(3);

maxBound_x = vmax(1);
maxBound_y = vmax(2);
maxBound_z = vmax(3);

    if (ray_direction_x >= 0) 
    	tmin = (minBound_x - ray_origin_x) / ray_direction_x;
    	tmax = (maxBound_x - ray_origin_x) / ray_direction_x;
    else
    	tmin = (maxBound_x - ray_origin_x) / ray_direction_x;
    	tmax = (minBound_x - ray_origin_x) / ray_direction_x;
    end
  
    if (ray_direction_y >= 0) 
        tymin = (minBound_y - ray_origin_y) / ray_direction_y;
        tymax = (maxBound_y - ray_origin_y) / ray_direction_y;
    else
    	tymin = (maxBound_y - ray_origin_y) / ray_direction_y;
    	tymax = (minBound_y - ray_origin_y) / ray_direction_y;
    end
    if ( (tmin > tymax) || (tymin > tmax) )
        flag = 0;
        tmin = -1;
    	return;
    end
       
    if (tymin > tmin)
        tmin = tymin;
    end
    
	if (tymax < tmax)
        tmax = tymax;
    end
    
	if (ray_direction_z >= 0)
       tzmin = (minBound_z - ray_origin_z) / ray_direction_z;
       tzmax = (maxBound_z - ray_origin_z) / ray_direction_z;
    else
       tzmin = (maxBound_z - ray_origin_z) / ray_direction_z;
       tzmax = (minBound_z - ray_origin_z) / ray_direction_z;
    end
    if ((tmin > tzmax) || (tzmin > tmax))
        flag = 0;
        tmin = -1;
       return;
    end
    
    if (tzmin > tmin)
        tmin = tzmin;
    end
   
    if (tzmax < tmax)
        tmax = tzmax;
    end
    
  % if( (tmin < t1) && (tmax > t0) )
      flag = 1;
  % else
  %    flag = 0;
  %    tmin = -1;
  % end;
end
