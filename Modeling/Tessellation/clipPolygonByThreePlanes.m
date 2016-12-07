function polygon_out = clipPolygonByThreePlanes(polygon_in, polygon_extent, bbox_half_width, plane_1, plane_2, plane_3)
polygon_out = polygon_in;
%disp('Three plane')

% Clip the polygon by the three planes
try
    polygon_out = clipConvexPolygon3dHP(clipConvexPolygon3dHP(clipConvexPolygon3dHP([planePoint(polygon_out, [-polygon_extent -polygon_extent]); planePoint(polygon_out, [-polygon_extent polygon_extent]); planePoint(polygon_out, [polygon_extent polygon_extent]); planePoint(polygon_out, [polygon_extent -polygon_extent])], plane_1), plane_2), plane_3);
catch
end
% Now clip by the bounding box
xp = createPlane([ bbox_half_width 0  0],  [1 0 0]);
xm = createPlane([-bbox_half_width 0  0], -[1 0 0]);
yp = createPlane([0  bbox_half_width  0],  [0 1 0]);
ym = createPlane([0 -bbox_half_width  0], -[0 1 0]);
zp = createPlane([0  0  bbox_half_width],  [0 0 1]);
zm = createPlane([0  0 -bbox_half_width], -[0 0 1]);
%
try
    %disp(size(clipConvexPolygon3dHP(polygon_out, xm)))
    if(size(clipConvexPolygon3dHP(polygon_out, xm), 1) ~= 0)
        try
            polygon_out = clipConvexPolygon3dHP(polygon_out, xm);
            %disp('Clipped')
        catch
            %disp('Error clipping')
        end
    end

    %disp(size(clipConvexPolygon3dHP(polygon_out, xp)))
    if(size(clipConvexPolygon3dHP(polygon_out, xp), 1) ~= 0)
        try
            polygon_out = clipConvexPolygon3dHP(polygon_out, xp);
            %disp('Clipped')
        catch
            %disp('Error clipping')
        end
    end

    %disp(size(clipConvexPolygon3dHP(polygon_out, ym)))
    if(size(clipConvexPolygon3dHP(polygon_out, yp), 1) ~= 0)
        try
            polygon_out = clipConvexPolygon3dHP(polygon_out, ym);
            %disp('Clipped')
        catch
            %disp('Error clipping')
        end
    end


    %disp(size(clipConvexPolygon3dHP(polygon_out, yp)))
    if(size(clipConvexPolygon3dHP(polygon_out, yp), 1) ~= 0)
        try
            polygon_out = clipConvexPolygon3dHP(polygon_out, yp);
            %disp('Clipped')
        catch
            %disp('Error clipping')
        end
    end


    %disp(size(clipConvexPolygon3dHP(polygon_out, zm)))
    if(size(clipConvexPolygon3dHP(polygon_out, zm), 1) ~= 0)
        try
            polygon_out = clipConvexPolygon3dHP(polygon_out, zm);
            %disp('Clipped')
        catch
            %disp('Error clipping')
        end
    end


    %disp(size(clipConvexPolygon3dHP(polygon_out, zp)))
    if(size(clipConvexPolygon3dHP(polygon_out, zp), 1) ~= 0)
        try
            polygon_out = clipConvexPolygon3dHP(polygon_out, zp);
            %disp('Clipped')
        catch
            %disp('Error clipping')
        end
    end
catch
end
% 
% try
%     polygon_out = clipConvexPolygon3dHP(polygon_out, xm);
%     disp('Clipped')
% catch
%     disp('Error clipping')
% end
% 
% try
%     polygon_out = clipConvexPolygon3dHP(polygon_out, xp);
%     disp('Clipped')
% catch
%     disp('Error clipping')
% end
% 
% try
%     polygon_out = clipConvexPolygon3dHP(polygon_out, ym);
%     disp('Clipped')
% catch
%     disp('Error clipping')
% end
% 
% try
%     polygon_out = clipConvexPolygon3dHP(polygon_out, yp);
%     disp('Clipped')
% catch
%     disp('Error clipping')
% end
% 
% try
%     polygon_out = clipConvexPolygon3dHP(polygon_out, zm);
%     disp('Clipped')
% catch
%     disp('Error clipping')
% end
% 
% try
%     polygon_out = clipConvexPolygon3dHP(polygon_out, zp);
%     disp('Clipped')
% catch
%     disp('Error clipping')
% end

end