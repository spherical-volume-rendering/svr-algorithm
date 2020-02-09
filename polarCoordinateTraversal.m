function polarCoordinateTraversal(origin, direction, polarGrid2D, verbose)
% Input:
%    ray origin.
%    ray direction.
%    polarGrid2D: grid dimensions (deltaRadius, deltaTheta, maxRadius).
%
% Pre-conditions: For simplcity, the ray's 'origin' is always along 'polarGrid2D.maxRadius'.
% Notes: CURRENTLY UNDER CONSTRUCTION.
        start   = origin; % + tmin*direction;
        boxSize = grid3D.maxBound-grid3D.minBound;
        
        if (verbose)
            plot3(start(1), start(2), start(3), 'r.', 'MarkerSize', 15);
        end;
        
        % Calculate the first voxel we are located in.
        x = floor( ((start(1)-grid3D.minBound(1))/boxSize(1))*grid3D.nx )+1;
        y = floor( ((start(2)-grid3D.minBound(2))/boxSize(2))*grid3D.ny )+1;
        z = floor( ((start(3)-grid3D.minBound(3))/boxSize(3))*grid3D.nz )+1;
        if (x==(grid3D.nx+1));  x=x-1;  end;
        if (y==(grid3D.ny+1));  y=y-1;  end;            
        if (z==(grid3D.nz+1));  z=z-1;  end;
        
        if (direction(1)>=0)
            tVoxelX = (x)/grid3D.nx;
            stepX = 1;
        else
            tVoxelX = (x-1)/grid3D.nx;
            stepX = -1;  
        end;
        
        if (direction(2)>=0)
            tVoxelY = (y)/grid3D.ny;
            stepY = 1;
        else
            tVoxelY = (y-1)/grid3D.ny;
            stepY = -1;
        end;
        
        if (direction(3)>=0)
            tVoxelZ = (z)/grid3D.nz; 
            stepZ = 1;
        else
            tVoxelZ = (z-1)/grid3D.nz;
            stepZ = -1;  
        end;
                
        voxelMaxX  = grid3D.minBound(1) + tVoxelX*boxSize(1);
        voxelMaxY  = grid3D.minBound(2) + tVoxelY*boxSize(2);
        voxelMaxZ  = grid3D.minBound(3) + tVoxelZ*boxSize(3);
        tMaxX      = tmin + (voxelMaxX-start(1))/direction(1);
        tMaxY      = tmin + (voxelMaxY-start(2))/direction(2);
        tMaxZ      = tmin + (voxelMaxZ-start(3))/direction(3);
        
        voxelSizeX = boxSize(1)/grid3D.nx;
        voxelSizeY = boxSize(2)/grid3D.ny;
        voxelSizeZ = boxSize(3)/grid3D.nz;        
        
        tDeltaX    = voxelSizeX/abs(direction(1));
        tDeltaY    = voxelSizeY/abs(direction(2));
        tDeltaZ    = voxelSizeZ/abs(direction(3));
                
        while ( (x<=grid3D.nx)&&(x>=1) && (y<=grid3D.ny)&&(y>=1) && (z<=grid3D.nz)&&(z>=1) )
            if (verbose)
                fprintf('\nIntersection: voxel = [%d %d %d]', [x y z]);
                
                t1 = [(x-1)/grid3D.nx, (y-1)/grid3D.ny, (z-1)/grid3D.nz ]';
                t2 = [  (x)/grid3D.nx,  (y)/grid3D.ny,    (z)/grid3D.nz ]';        
                vmin = (grid3D.minBound + t1.*boxSize)';
                vmax = (grid3D.minBound + t2.*boxSize)';
                smallBoxVertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
                smallBoxFaces    = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
 
                h = patch('Vertices', smallBoxVertices, 'Faces', smallBoxFaces, 'FaceColor', 'blue', 'EdgeColor', 'white');
                set(h,'FaceAlpha',0.2);
            end;
            
            % ---------------------------------------------------------- %
            % check if voxel [x,y,z] contains any intersection with the ray
            %
            %   if ( intersection )
            %       break;
            %   end;
            % ---------------------------------------------------------- %
            
            if (tMaxX < tMaxY)
                if (tMaxX < tMaxZ)
                    x = x + stepX;
                    tMaxX = tMaxX + tDeltaX;
                else
                    z = z + stepZ;
                    tMaxZ = tMaxZ + tDeltaZ;
                end;
            else
                if (tMaxY < tMaxZ)
                    y = y + stepY;
                    tMaxY = tMaxY + tDeltaY;             
                else
                    z = z + stepZ;
                    tMaxZ = tMaxZ + tDeltaZ;
                end;
            end;
        end;        
     end;
end
