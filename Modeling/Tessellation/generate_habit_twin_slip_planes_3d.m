function [fighandle, elements, element_phase] = generate_habit_twin_slip_planes_3d(austenite_color_ip, CV1_color_ip, CV2_color_ip, ...
                                                      habit_plane_normal_ip, twin_plane_normal_ip, ...
                                                      CV_volume_ratio_ip, slip_plane_normal, slip_direction, ...
                                                      slip_plane_color_ip, slip_plane_opacity, foil_plane, ...
                                                      bbox_half_width, number_of_twins, ...
                                                      mesh_data_dir, node_file_name, connectivity_file_name, ...
                                                      elems_input, plot_flag)
% 3D plot showing habit and twin interafces for CV 1-9.
% CHESS June 2015 SC data
%
% Specify colors for the phases
austenite_color    = austenite_color_ip;
CV1_color          = CV1_color_ip;
CV2_color          = CV2_color_ip;
habit_plane_normal = habit_plane_normal_ip;
twin_plane_normal  = twin_plane_normal_ip;
CV_volume_ratio    = CV_volume_ratio_ip;
slip_direction     = slip_direction / norm(slip_direction);
%
fighandle = figure('Position', [200 200 600 600]);
% Hack to force 3D view. Not sure if needed in 2015a.
if(plot_flag > 0)
    plot3(0, 0, 0);
end
% Various planes are defined by their normal. a variable ending in 'p'
% means that the normal is +ve. 'm' at the end means that the normal is -ve
% So planem and planep define the same plane but with the normal pointing
% in the opposite direction.
%
% Foil plane
if(~isempty(foil_plane))
    foil = createPlane([0 0 0], foil_plane);
end
% Invariant/habit plane
i = createPlane([0 0 0], -habit_plane_normal);
im = createPlane([0 0 0], habit_plane_normal);
% Twin plane
t = createPlane([0 0 0], twin_plane_normal);
tm = createPlane([0 0 0], -twin_plane_normal);
% It's neat to show multiple twin planes separated by the right distance.
% This distance is taken as proportional to the volume fraction ratio of
% CV1 and CV2.
CV_volume_ratio_prime = 1.0 - CV_volume_ratio;
twin_offset = 0;
% Larger multiplier means fewer twin pairs will be plotted.
twin_spacing_multiplier = 5/number_of_twins;
%
t_store = cell(4 * number_of_twins + 1, 1);
tm_store = cell(4 * number_of_twins + 1, 1);
t_store_tmp = cell(4 * number_of_twins + 1, 1);
tm_store_tmp = cell(4 * number_of_twins + 1, 1);
t_offset = 0;
t_offset_store = zeros(4 * number_of_twins + 1, 1);
for ii = 1:number_of_twins
    t_store_tmp{ii}                          = parallelPlane(t, twin_offset + twin_spacing_multiplier*(t_offset + CV_volume_ratio));
    t_offset_store(ii)                       = t_offset + CV_volume_ratio;
    t_store_tmp{1 * number_of_twins + ii}    = parallelPlane(t, twin_offset + twin_spacing_multiplier*(t_offset + 1));
    t_offset_store(1 * number_of_twins + ii) = t_offset + 1;
    t_store_tmp{2 * number_of_twins + ii}    = parallelPlane(t, twin_offset - twin_spacing_multiplier*(t_offset + 1));
    t_offset_store(2 * number_of_twins + ii) = -t_offset - 1;
    t_store_tmp{3 * number_of_twins + ii}    = parallelPlane(t, twin_offset - twin_spacing_multiplier*(t_offset + CV_volume_ratio_prime));
    t_offset_store(3 * number_of_twins + ii) = -t_offset - CV_volume_ratio_prime;
    tm_store_tmp{ii}                         = parallelPlane(tm, -twin_offset - twin_spacing_multiplier*(t_offset + CV_volume_ratio));
    tm_store_tmp{1 * number_of_twins + ii}   = parallelPlane(tm, -twin_offset - twin_spacing_multiplier*(t_offset + 1));
    tm_store_tmp{2 * number_of_twins + ii}   = parallelPlane(tm, -twin_offset + twin_spacing_multiplier*(t_offset + 1));
    tm_store_tmp{3 * number_of_twins + ii}   = parallelPlane(tm, -twin_offset + twin_spacing_multiplier*(t_offset + CV_volume_ratio_prime));
    t_offset = t_offset + 1;
end
t_store_tmp{end} = t;
tm_store_tmp{end} = tm;
t_offset_store(end) = 0.0;
% Sort the planes according to their offset (easier for plotting later)
[~, idx] = sort(t_offset_store);
for ii = 1:(4 * number_of_twins + 1)
    t_store{ii} = t_store_tmp{idx(ii)};
    tm_store{ii} = tm_store_tmp{idx(ii)};
end

% Slip plane
if(~isempty(slip_plane_normal) || ~isempty(slip_direction))
    s1 = createPlane([0 0 0],  slip_plane_normal);
    s2 = createPlane([0.35 0.35 0.35],  slip_plane_normal);
    s3 = createPlane([-0.35 -0.35 -0.35],  slip_plane_normal);
    s4 = createPlane([0.7 0.7 0.7],  slip_plane_normal);
    s5 = createPlane([-0.7 -0.7 -0.7],  slip_plane_normal);
    s1m = createPlane([0 0 0],  -slip_plane_normal);
    s2m = createPlane([0.35 0.35 0.35],  -slip_plane_normal);
    s3m = createPlane([-0.35 -0.35 -0.35],  -slip_plane_normal);
end
% Planes along X, Y, Z faces
xp = createPlane([ bbox_half_width 0  0], -[1 0 0]);
xm = createPlane([-bbox_half_width 0  0],  [1 0 0]);
yp = createPlane([0  bbox_half_width  0], -[0 1 0]);
ym = createPlane([0 -bbox_half_width  0],  [0 1 0]);
zp = createPlane([0  0  bbox_half_width], -[0 0 1]);
zm = createPlane([0  0 -bbox_half_width],  [0 0 1]);
% Plot stuff
% Invariant plane  = green
% Twin planes = magenta
if(plot_flag > 0)
    hold on;
end
% plane_extent is a terrible hack. plane_extent specifies the "width" of a
% plane. This is set to a very large number compared to our plot bounding
% box. That way, any two reasonably non-parallel parallel planes are
% guaranteed to intersect and geom3d does not throw an exception. This hack
% pretty much prevents exporting the model to VRML, U3D etc, unless we come
% up with an easy way of trimming the planes outside the bounding box
% (2x2x2 for now)
plane_extent = 1800;
% Foil plane
if(~isempty(foil_plane) && plot_flag > 0)
    try
        patch_foil = fillPolygon3d(clipPolygonByBBox(foil, plane_extent, bbox_half_width), 'w');
    catch
    end
end
% Habit plane
if(plot_flag > 0)
    patch_a = fillPolygon3d(clipPolygonByBBox(i, plane_extent, bbox_half_width), austenite_color);
end
% Clipped bounding box for the habit plane
if((isempty(slip_plane_normal) || isempty(slip_direction)) && plot_flag > 0)
    try
        fillPolygon3d(clipPolygonByOnePlanes(xp, plane_extent, bbox_half_width, im), austenite_color);
        fillPolygon3d(clipPolygonByOnePlanes(xm, plane_extent, bbox_half_width, im), austenite_color);
        fillPolygon3d(clipPolygonByOnePlanes(yp, plane_extent, bbox_half_width, im), austenite_color);
        fillPolygon3d(clipPolygonByOnePlanes(ym, plane_extent, bbox_half_width, im), austenite_color);
        fillPolygon3d(clipPolygonByOnePlanes(zp, plane_extent, bbox_half_width, im), austenite_color);
        fillPolygon3d(clipPolygonByOnePlanes(zm, plane_extent, bbox_half_width, im), austenite_color);
    catch
    end
end
% Clipped slip planes
if((~isempty(slip_plane_normal) || ~isempty(slip_direction)) && plot_flag > 0)
    try
        patch_slip_1 = fillPolygon3d(clipPolygonByOnePlanes(s1, plane_extent, bbox_half_width, im), slip_plane_color_ip);
        patch_slip_2 = fillPolygon3d(clipPolygonByOnePlanes(s2, plane_extent, bbox_half_width, im), slip_plane_color_ip);
        patch_slip_3 = fillPolygon3d(clipPolygonByOnePlanes(s3, plane_extent, bbox_half_width, im), slip_plane_color_ip);
        patch_slip_4 = fillPolygon3d(clipPolygonByOnePlanes(s4, plane_extent, bbox_half_width, im), slip_plane_color_ip);
        patch_slip_5 = fillPolygon3d(clipPolygonByOnePlanes(s5, plane_extent, bbox_half_width, im), slip_plane_color_ip);
        set(patch_slip_1, 'FaceAlpha', slip_plane_opacity);
        set(patch_slip_2, 'FaceAlpha', slip_plane_opacity);
        set(patch_slip_3, 'FaceAlpha', slip_plane_opacity);
        set(patch_slip_4, 'FaceAlpha', slip_plane_opacity);
        set(patch_slip_5, 'FaceAlpha', slip_plane_opacity);
    catch
    end
end
% Clipped twin planes
t_drawn = zeros(4 * number_of_twins + 2, 1);
for ii = 1:(4 * number_of_twins + 1)
    try
        if(plot_flag > 0)
            fillPolygon3d(clipPolygonByOnePlanes(t_store{ii}, plane_extent, bbox_half_width, i), 'm');
            t_drawn(ii) = 1;
        end
    catch
        % disp(['Did not plot twin plane ' num2str(ii)])
    end
end

% Clipped bounding box for the twins
%
% +ve X face
for ii = 2:2:(4 * number_of_twins)
    try
        if(plot_flag > 0)
            patch_cv1 = fillPolygon3d(clipPolygonByThreePlanes(xp, plane_extent, bbox_half_width, i, tm_store{ii - 1}, t_store{ii}), CV1_color);
        end
    catch
    end
    try
        if(plot_flag > 0)
            fillPolygon3d(clipPolygonByThreePlanes(xp, plane_extent, bbox_half_width, i, tm_store{ii}, t_store{ii + 1}), CV2_color);
        end
    catch
    end
end
%
% -ve X face
for ii = 2:2:(4 * number_of_twins)
    try
        if(plot_flag > 0)
            patch_cv1 = fillPolygon3d(clipPolygonByThreePlanes(xm, plane_extent, bbox_half_width, i, tm_store{ii - 1}, t_store{ii}), CV1_color);
        end
    catch
    end
    try
        if(plot_flag > 0)
            fillPolygon3d(clipPolygonByThreePlanes(xm, plane_extent, bbox_half_width, i, tm_store{ii}, t_store{ii + 1}), CV2_color);
        end
    catch
    end
end
%
% +ve Y face
for ii = 2:2:(4 * number_of_twins)
    try
        if(plot_flag > 0)
            patch_cv1 = fillPolygon3d(clipPolygonByThreePlanes(yp, plane_extent, bbox_half_width, i, tm_store{ii - 1}, t_store{ii}), CV1_color);
        end
    catch
    end
    try
        if(plot_flag > 0)
            fillPolygon3d(clipPolygonByThreePlanes(yp, plane_extent, bbox_half_width, i, tm_store{ii}, t_store{ii + 1}), CV2_color);
        end
    catch
    end
end
%
% -ve Y face
for ii = 2:2:(4 * number_of_twins)
    try
        if(plot_flag > 0)
            patch_cv1 = fillPolygon3d(clipPolygonByThreePlanes(ym, plane_extent, bbox_half_width, i, tm_store{ii - 1}, t_store{ii}), CV1_color);
        end
    catch
    end
    try
        if(plot_flag > 0)
            fillPolygon3d(clipPolygonByThreePlanes(ym, plane_extent, bbox_half_width, i, tm_store{ii}, t_store{ii + 1}), CV2_color);
        end
    catch
    end
end
%
% +ve Z face
for ii = 2:2:(4 * number_of_twins)
    try
        if(plot_flag > 0)
            patch_cv1 = fillPolygon3d(clipPolygonByThreePlanes(zp, plane_extent, bbox_half_width, i, tm_store{ii - 1}, t_store{ii}), CV1_color);
        end
    catch
    end
    try
        if(plot_flag > 0)
            fillPolygon3d(clipPolygonByThreePlanes(zp, plane_extent, bbox_half_width, i, tm_store{ii}, t_store{ii + 1}), CV2_color);
        end
    catch
    end
end
%
% -ve Z face
for ii = 2:2:(4 * number_of_twins)
    try
        if(plot_flag > 0)
            patch_cv1 = fillPolygon3d(clipPolygonByThreePlanes(zm, plane_extent, bbox_half_width, i, tm_store{ii - 1}, t_store{ii}), CV1_color);
        end
    catch
    end
    try
        if(plot_flag > 0)
            fillPolygon3d(clipPolygonByThreePlanes(zm, plane_extent, bbox_half_width, i, tm_store{ii}, t_store{ii + 1}), CV2_color);
        end
    catch
    end
end
%
% Slip direction
% line([0 0 + slip_direction(1)], [1.5 1.5 + slip_direction(2)], [0 0 + slip_direction(3)], 'LineWidth', 3, 'Color', 'k');
% X, Y, Z axes
if(plot_flag > 0)
    line([-2 -1],   [-2 -2],   [-2 -2],   'Color', 'r', 'LineWidth', 3)
    line([-2 -2],   [-2 -1],   [-2 -2],   'Color', 'g', 'LineWidth', 3)
    line([-2 -2],   [-2 -2],   [-2 -1],   'Color', 'b', 'LineWidth', 3)
    hold off;
end
% Set axis etc
if(plot_flag > 0)
    alpha(0.8);
    % xlim([-bbox_half_width-1 bbox_half_width+1]); ylim([-bbox_half_width-1 bbox_half_width+1]); zlim([-bbox_half_width-1 bbox_half_width+1]);
    xlim([-bbox_half_width bbox_half_width]); ylim([-bbox_half_width bbox_half_width]); zlim([-bbox_half_width bbox_half_width]);
    grid off; box on; axis vis3d;
    set(gca, 'XTick', []); set(gca, 'XTickLabel', []);
    set(gca, 'YTick', []); set(gca, 'YTickLabel', []);
    set(gca, 'ZTick', []); set(gca, 'ZTickLabel', []);
    camproj('perspective')
    camtarget([-1 1 1]);
    camup([0 1 0]);
    campos([10 30 30]);
    % Lighting
    camlight;
    set(findall(gcf, 'Type', 'patch'), 'DiffuseStrength', 0.8)
    set(findall(gcf, 'Type', 'patch'), 'AmbientStrength', 0.6)
    % Legend
    try
        l = legend([patch_a, patch_cv1, patch_cv2], 'Austenite', 'CV 1', 'CV 2');
        set(l, 'FontSize', 14);
    catch
    end
    %
    % Some of the clipped patches are still going to be out of the bounding
    % box. Here we prune those.
    p = findall(gcf, 'Type', 'patch');
    for ii = 1:size(p(:))
        p_centroid = mean(p(ii).Vertices, 1);
        if(p_centroid(1) > bbox_half_width || p_centroid(2) > bbox_half_width || ...
                p_centroid(3) > bbox_half_width || p_centroid(1) < -bbox_half_width || ...
                p_centroid(2) < -bbox_half_width || p_centroid(3) < -bbox_half_width)
            delete(p(ii));
        end
    end
    %
    if(~isempty(foil_plane))
        set(patch_foil, 'FaceAlpha', 0.3);
    end
end
%
% Create a 3D grid and assign appropriate phase to the point (1 = CV1, 2 = Austenite, 3 = CV 2)
% If the nodes/connectivity file names are empty, create a uniform grid
if((isempty(mesh_data_dir) || isempty(node_file_name) || isempty(connectivity_file_name)) && isempty(elems_input))
    % Uniform grid
    [X, Y, Z] = ndgrid(-bbox_half_width:0.1:bbox_half_width, -bbox_half_width:0.1:bbox_half_width, -bbox_half_width:0.1:bbox_half_width);
    elements = [X(:) Y(:) Z(:)];
    user_grid_bbox_ratio = 1.0;
    user_grid_centroid = [0 0 0];
elseif(~isempty(elems_input))
    % User supplied grid in a variable (Mx3: coordinates for M elements)
    elements = elems_input;
    user_grid_centroid = mean(elements, 1);
    elements = elements - repmat(user_grid_centroid, size(elements, 1), 1);
    user_grid_bbox_ratio = max(elements(:)) / bbox_half_width;
    elements = elements / user_grid_bbox_ratio;
else
    % User supplied grid in Abaqus nodes and connectivity format
    nodes_in = importdata(fullfile(mesh_data_dir, node_file_name));
    connectivity_in = importdata(fullfile(mesh_data_dir, connectivity_file_name));
    nodes(nodes_in(:, 1), :) = nodes_in(:, 2:end);
    connectivity(connectivity_in(:, 1), :) = connectivity_in(:, 2:end);
    elements = zeros(size(connectivity, 1), 3);
    for ii = 1:size(elements, 1)
        elements(ii, :) = mean(nodes(connectivity(ii, :)', :), 1);
    end
    % User  supplied grid can have any bbox. However we are using a
    % different and symmetric bbox for the visualizaion. Since the 
    % microstructure is made using the visualization, we need to 
    % temporarily make the bboxes equal. We will return elems corresponding
    % to whatever grid data the user supplied.
    user_grid_centroid = mean(elements, 1);
    elements = elements - repmat(user_grid_centroid, size(elements, 1), 1);
    user_grid_bbox_ratio = max(elements(:)) / bbox_half_width;
    elements = elements / user_grid_bbox_ratio;
end
element_phase = 3 * ones(size(elements, 1), 1);
% Determine which points belong to a CV
for ii = 1:(4 * number_of_twins)
    element_phase(~isBelowPlane(elements, t_store{ii + 1}) & isBelowPlane(elements, tm_store{ii})) = ii;
end
element_phase(mod(element_phase, 2) == 1) = 1;
element_phase(mod(element_phase, 2) == 0) = 3;
% Which points are below the invariant plane
is_below_i = isBelowPlane(elements, i);
element_phase(~is_below_i) = 2;
%
% Plot
if(plot_flag > 0)
    figure('Position', [100, 100, 600, 600]);
    scatter3(elements(:, 1), elements(:, 2), elements(:, 3), 36, element_phase, 'filled');
    colormap jet;
    axis equal; axis vis3d;
    camproj('perspective')
    camtarget([-1 1 1]);
    camup([0 1 0]);
    campos([10 30 30]);
    % All done. Now change grid bbox to what te user supplied.
    elements = elements * user_grid_bbox_ratio;
    elements = elements + repmat(user_grid_centroid, size(elements, 1), 1);
end
end
