function polygon_out = clipPolygonByBBox(polygon_in, polygon_extent, bbox_half_width)
polygon_out = polygon_in;

%disp('BBox')

polygon_out = [planePoint(polygon_out, [-polygon_extent -polygon_extent]); planePoint(polygon_out, [-polygon_extent polygon_extent]); planePoint(polygon_out, [polygon_extent polygon_extent]); planePoint(polygon_out, [polygon_extent -polygon_extent])];

% Clip by the bounding box
xp = createPlane([ bbox_half_width 0  0],  [1 0 0]);
xm = createPlane([-bbox_half_width 0  0], -[1 0 0]);
yp = createPlane([0  bbox_half_width  0],  [0 1 0]);
ym = createPlane([0 -bbox_half_width  0], -[0 1 0]);
zp = createPlane([0  0  bbox_half_width],  [0 0 1]);
zm = createPlane([0  0 -bbox_half_width], -[0 0 1]);
%
try
    polygon_out = clipConvexPolygon3dHP(polygon_out, xm);
    %disp('Clipped')
catch
    %disp('Error clipping')
end

try
    polygon_out = clipConvexPolygon3dHP(polygon_out, xp);
    %disp('Clipped')
catch
    %disp('Error clipping')
end

try
    polygon_out = clipConvexPolygon3dHP(polygon_out, ym);
    %disp('Clipped')
catch
    %disp('Error clipping')
end

try
    polygon_out = clipConvexPolygon3dHP(polygon_out, yp);
    %disp('Clipped')
catch
    %disp('Error clipping')
end

try
    polygon_out = clipConvexPolygon3dHP(polygon_out, zm);
    %disp('Clipped')
catch
    %disp('Error clipping')
end

try
    polygon_out = clipConvexPolygon3dHP(polygon_out, zp);
    %disp('Clipped')
catch
    %disp('Error clipping')
end

end