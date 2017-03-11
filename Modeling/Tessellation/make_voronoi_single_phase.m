function [tess_coords, tess, tess_quat, strain, intensity_multiplier] = make_voronoi_single_phase(generate_random_grains, ...
                                                                        num_of_grains, sigma_o, mag_e, sigma_e, grains_com, ...
                                                                        quat, strain, ...
                                                                        read_grid_data, ...
                                                                        sample_size_um, gridpoints_per_dim, ...
                                                                        gridpoints_coord, ...
                                                                        ms_file_name, ...
                                                                        add_powder_to_gb, powder_intensity_multiplier, ...
                                                                        material_name, ...
                                                                        add_local_perturbations, simulate_intensity_falloff_at_gb)
%% Synthetic microstructure generation: Single phase
%
% PREREQUISITES: MTEX and export_fig packages are required.
% INPUTS:
% generate_random_grains      = 1: generate random grain COM. Must specify grain number, sigma_o, mag_e, sigma_e.
% num_of_grains               = Number of grains (required for both random grains and user-specified grain COMs)
% sigma_o                     = Amount of orientation spread (mosaicity) to add to each grain
% mag_e                       = Magnitude of the strain in the grain. (A random number between 0 and mag_e will be assigned to each strain component)
% sigma_e                     = Amount of deformation spread to add to each grain
% grains_com                  = Array of (num_of_grains, 3) size
%                               specifying the coordinates of grain centroid.
% quat                        = Orientation of each grain in quaternion
%                               format. Size (num_of_grains, 4)
% strain                      = strain tensor in each grain. Size
%                               (num_of_grains, 6) [components: 11, 22, 33, 12, 13, 23]
% read_grain_data             = 1: read mesh grid points from the matrix
%                               gridpoints_coord. 0: Generate a random uniform grid with density
%                               gridpoints_per_dim along each direction
% ms_file_name                = Name of the file to which the output is
%                               written. This file will be the input for forward modeling code. 
%                               Default: 'ms-synth-data.csv'
%%
if(generate_random_grains == 0)
    % DO NOT generate random grains
    % Get grain number and grain centroid data from the user
    com = grains_com;                        % Get grain center of mass data from the user.
    M   = size(com, 1);                      % Get number of grains
    assert(size(com, 1) == size(quat, 1));   % A quaternion (orientation) must be specified for each grain
    assert(size(com, 1) == size(strain, 1)); % A strain array (6 elements) must be specified for each grain
else
    % Get the number of random grains to generate from the user.
    % Generate random COM, orientation, strain data for each grain.
    L      = sample_size_um;          % Sample size in microns
    M      = num_of_grains;           % Get number of grains from user
    com    = (rand(M, 3) - 0.5) * L;  % Random grain centroid coordinates
    quat   = randq(M);                % Random quaternions for grain orientations
    quat   = squeeze(double(quat));   % Convert from quat MTEX object to normal matrix
    strain = rand(M, 6) * mag_e;      % Small random strain at each grain
end
L = sample_size_um;                   % Sample dimensions in microns
if(read_grid_data == 0)
    % DO NOT read grid data. Create a uniform grid with this density.
    N = gridpoints_per_dim;           % NxNxN uniform grid will be created from the specimen
else
    N = floor(1.5 * power(size(gridpoints_coord, 1), 1.0/3.0));
end
%
% Create a tessellation
tess = zeros(N, N, N, 1);             % A 3D matrix containing grain ID for each point
tess_coords = zeros(N*N*N, 3);        % A 2D matrix containing coordinates of each point
tess_quat = zeros(N, N, N, 4);        % A 4D matrix containing four quaternion components at each tessellation point
for ii = 1:N
    for jj = 1:N
        for kk = 1:N
            tess_coords(N*N*(ii - 1) + N*(jj - 1) + kk, :) = [ii/N*L - L/2.0 jj/N*L - L/2.0 kk/N*L - L/2.0];
        end
    end
end
%
% If output microstructure filename is empty, set a default.
if(isempty(ms_file_name))
    ms_file_name = 'ms-synth-data.csv';
end
% If no name is specified for the material, set a default.
if(isempty(material_name))
    material_name = 'NiTi_cubic';
end
%
% Intensity multiplier value for the powder data at added at grain boundaries
if(isempty(powder_intensity_multiplier) && add_powder_to_gb == 1)
    powder_intensity_multiplier = 0.1;
end
%%
% If single crystal, make sure quat is still a row vector.
if(M == 1)
    quat = reshape(quat, [1 4]);
end
%
% Calculate distance between each material point and grain centroids
dist = pdist2(tess_coords, com);
%
for ii = 1:N
    for jj = 1:N
        for kk = 1:N
            % Voronoi tessellation is where each material point is assigned
            % to a grain whose COM is closest to that material point
            [~, tess(ii, jj, kk)] = min(dist(N*N*(ii - 1) + N*(jj - 1) + kk, :));
            tess_quat(ii, jj, kk, 1) = quat(tess(ii, jj, kk), 1);
            tess_quat(ii, jj, kk, 2) = quat(tess(ii, jj, kk), 2);
            tess_quat(ii, jj, kk, 3) = quat(tess(ii, jj, kk), 3);
            tess_quat(ii, jj, kk, 4) = quat(tess(ii, jj, kk), 4);
        end
    end
end
%%
% Calculate intensity multiplier. This is a hack to simulate the diffuse
% nature of spots. We define a bounding box of fixed size (approx. N/10).
% Move the box to a parent material point. Calculate the number of points
% inside that box that have the same orientation as the parent. Normalize
% this number by the box volume. This procedure will give a large number
% for a parent point that is near the grain center. Near the grain
% boundaries, the number will drop significantly. Forward modeling program
% will use this number, multiply the fixed intensity with it and store the
% resultant intensity value in the GE2 file. Thus points near grain
% boundaries will make smaller contribution compared to the points near
% grain center.
intensity_multiplier = zeros(N, N, N);
bbox_size = round([N N N]/4);
for ii = 1:3
    if(mod(bbox_size(ii), 2) == 1)
        bbox_size(ii) = bbox_size(ii) + 1;
    end
end
if(simulate_intensity_falloff_at_gb ~= 0)
    parfor ii = 1:N
        for jj = 1:N
            for kk = 1:N
                bbox_min_1 = max(ii - bbox_size(1)/2, 1);
                bbox_max_1 = min(ii + bbox_size(1)/2, N);
                bbox_min_2 = max(jj - bbox_size(2)/2, 1);
                bbox_max_2 = min(jj + bbox_size(2)/2, N);
                bbox_min_3 = max(kk - bbox_size(3)/2, 1);
                bbox_max_3 = min(kk + bbox_size(3)/2, N);
                %
                bbox_grain_ids = tess(bbox_min_1:bbox_max_1, ...
                    bbox_min_2:bbox_max_2, ...
                    bbox_min_3:bbox_max_3);
                bbox_same_id_count = sum(bbox_grain_ids(:) == tess(ii, jj, kk));
                intensity_multiplier(ii, jj, kk) = bbox_same_id_count / (bbox_size(1) * bbox_size(2) * bbox_size(3));
            end
        end
    end
    intensity_multiplier = intensity_multiplier / max(intensity_multiplier(:));
else
    intensity_multiplier = ones(N, N, N);
end
%%
local_perturbations = rand(N, N, N, 10);
quat_perturbation_norm = zeros(N, N, N);
stretch_perturbation_norm = zeros(N, N, N);
%
local_perturbations_tmp = local_perturbations;
parfor ii = 1:N
    for jj = 1:N
        for kk = 1:N
            bbox_min_1 = max(ii - bbox_size(1)/2, 1);
            bbox_max_1 = min(ii + bbox_size(1)/2, N);
            bbox_min_2 = max(jj - bbox_size(2)/2, 1);
            bbox_max_2 = min(jj + bbox_size(2)/2, N);
            bbox_min_3 = max(kk - bbox_size(3)/2, 1);
            bbox_max_3 = min(kk + bbox_size(3)/2, N);
            %
            bbox_perturb = local_perturbations_tmp(bbox_min_1:bbox_max_1, ...
                bbox_min_2:bbox_max_2, ...
                bbox_min_3:bbox_max_3);
            local_perturbations(ii, jj, kk) = mean(bbox_perturb(:));
        end
    end
end
% Add powder at grain boundaries
if(add_powder_to_gb)
    bbox_size = round([N N N]/10);
    quat_powder   = randq(10000);                   % Random quaternions for powder data
    quat_powder   = squeeze(double(quat_powder));   % Convert from quat MTEX object to normal matrix
    is_gb         = zeros(N, N, N);
    for ii = 1:3
        if(mod(bbox_size(ii), 2) == 1)
            bbox_size(ii) = bbox_size(ii) + 1;
        end
    end
    parfor ii = 1:N
        for jj = 1:N
            for kk = 1:N
                bbox_min_1 = max(ii - bbox_size(1)/2, 1);
                bbox_max_1 = min(ii + bbox_size(1)/2, N);
                bbox_min_2 = max(jj - bbox_size(2)/2, 1);
                bbox_max_2 = min(jj + bbox_size(2)/2, N);
                bbox_min_3 = max(kk - bbox_size(3)/2, 1);
                bbox_max_3 = min(kk + bbox_size(3)/2, N);
                %
                bbox_grain_ids = tess(bbox_min_1:bbox_max_1, ...
                    bbox_min_2:bbox_max_2, ...
                    bbox_min_3:bbox_max_3);
                bbox_unique_ids = unique(bbox_grain_ids(:));
                is_this_at_gb = 1;
                for mm = 1:numel(bbox_unique_ids)
                    if(sum(bbox_grain_ids(:) == bbox_unique_ids(mm)) / (bbox_size(1) * bbox_size(2) * bbox_size(3)) > 0.75)
                        is_this_at_gb = 0;
                    end
                end
                %
                if(is_this_at_gb == 1)
                    tess_quat(ii, jj, kk, :) = quat_powder(randi(10000), :);
                    intensity_multiplier(ii, jj, kk) = powder_intensity_multiplier;
                    is_gb(ii, jj, kk) = 1;
                end
            end
        end
    end
else
    is_gb         = zeros(N, N, N);
end
%%
% Generate output data at appropriate grid points. If the user did not
% speicify gridpoints_coord, then we just output the data at the uniform
% grid points that were created. Otherwise, we have to interpolate quat,
% strain and intensity_multiplier to gridpoints_coord using the
% nearest-neighbor interpolation scheme.
if(isempty(gridpoints_coord))
    % File to write microstructural data
    f = fopen(ms_file_name, 'w');
    for ii = 1:N
        for jj = 1:N
            for kk = 1:N
                % Write microstructure data in the column format that fwdmodel
                % python code reads.
                % Strain is really right stretch tensor.
                % I am adding a small noise to orientations and strains. Hence
                % the rand() calls.
                fprintf(f, '%12.6f, %12.6f, %12.6f, %s, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f\n', ...
                    ii/N*L - L/2.0, ...
                    jj/N*L - L/2.0, ...
                    kk/N*L - L/2.0, ...
                    material_name, ...
                    tess_quat(ii, jj, kk, 1) + local_perturbations(ii, jj, kk, 1)*sigma_o, ...
                    tess_quat(ii, jj, kk, 2) + local_perturbations(ii, jj, kk, 2)*sigma_o, ...
                    tess_quat(ii, jj, kk, 3) + local_perturbations(ii, jj, kk, 3)*sigma_o, ...
                    tess_quat(ii, jj, kk, 4) + local_perturbations(ii, jj, kk, 4)*sigma_o, ...
                    1.0 + strain(tess(ii, jj, kk), 1) + local_perturbations(ii, jj, kk, 5)*sigma_e, ...
                    1.0 + strain(tess(ii, jj, kk), 2) + local_perturbations(ii, jj, kk, 6)*sigma_e, ...
                    1.0 + strain(tess(ii, jj, kk), 3) + local_perturbations(ii, jj, kk, 7)*sigma_e, ...
                    strain(tess(ii, jj, kk), 4) + local_perturbations(ii, jj, kk, 8)*sigma_e, ...
                    strain(tess(ii, jj, kk), 5) + local_perturbations(ii, jj, kk, 9)*sigma_e, ...
                    strain(tess(ii, jj, kk), 6) + local_perturbations(ii, jj, kk, 10)*sigma_e, ...
                    intensity_multiplier(ii, jj, kk));
            end
        end
    end
    % Don't forget to close files.
    fclose(f);
else
    % Do nearest-neighbor interpolation
    output_tess                      = zeros(size(gridpoints_coord, 1), 1);
    output_tess_quat_1               = zeros(size(gridpoints_coord, 1), 1);
    output_tess_quat_2               = zeros(size(gridpoints_coord, 1), 1);
    output_tess_quat_3               = zeros(size(gridpoints_coord, 1), 1);
    output_tess_quat_4               = zeros(size(gridpoints_coord, 1), 1);
    output_tess_intensity_multiplier = zeros(size(gridpoints_coord, 1), 1);
    %
    tess_quat_reshaped = reshape(tess_quat, N*N*N, 4);
    %
    xgrid_vector = (-L/2.0):(L/(N - 1)):(L/2.0);
    ygrid_vector = (-L/2.0):(L/(N - 1)):(L/2.0);
    zgrid_vector = (-L/2.0):(L/(N - 1)):(L/2.0);
    [xvector_grid, yvector_grid, zvector_grid] = meshgrid(xgrid_vector, ygrid_vector, zgrid_vector);
    %
    parfor ii = 1:size(gridpoints_coord, 1)
        output_tess(ii) = interp3(xvector_grid, yvector_grid, zvector_grid, ...
                                  tess, gridpoints_coord(ii, 1), gridpoints_coord(ii, 2), gridpoints_coord(ii, 3), ...
                                  'nearest');
%         output_tess_quat_1(ii)    = interp3(xvector_grid, yvector_grid, zvector_grid, ...
%                                   squeeze(tess_quat(:, :, :, 1)), gridpoints_coord(ii, 1), gridpoints_coord(ii, 2), gridpoints_coord(ii, 3), ...
%                                   'nearest');
%         output_tess_quat_2(ii)    = interp3(xvector_grid, yvector_grid, zvector_grid, ...
%                                   squeeze(tess_quat(:, :, :, 2)), gridpoints_coord(ii, 1), gridpoints_coord(ii, 2), gridpoints_coord(ii, 3), ...
%                                   'nearest');
%         output_tess_quat_3(ii)    = interp3(xvector_grid, yvector_grid, zvector_grid, ...
%                                   squeeze(tess_quat(:, :, :, 3)), gridpoints_coord(ii, 1), gridpoints_coord(ii, 2), gridpoints_coord(ii, 3), ...
%                                   'nearest');
%         output_tess_quat_4(ii)    = interp3(xvector_grid, yvector_grid, zvector_grid, ...
%                                   squeeze(tess_quat(:, :, :, 4)), gridpoints_coord(ii, 1), gridpoints_coord(ii, 2), gridpoints_coord(ii, 3), ...
%                                   'nearest');
        output_tess_quat_1(ii)    = quat(output_tess(ii), 1);
        output_tess_quat_2(ii)    = quat(output_tess(ii), 2);
        output_tess_quat_3(ii)    = quat(output_tess(ii), 3);
        output_tess_quat_4(ii)    = quat(output_tess(ii), 4);
        if(simulate_intensity_falloff_at_gb ~= 0)
            output_tess_intensity_multiplier(ii) = interp3(xvector_grid, yvector_grid, zvector_grid, ...
                                      intensity_multiplier, gridpoints_coord(ii, 1), gridpoints_coord(ii, 2), gridpoints_coord(ii, 3), ...
                                      'nearest');
        else
            output_tess_intensity_multiplier(ii) = 1.0
        end
        %
        if(add_powder_to_gb)
            output_is_gb(ii) = interp3(xvector_grid, yvector_grid, zvector_grid, ...
                                      is_gb, gridpoints_coord(ii, 1), gridpoints_coord(ii, 2), gridpoints_coord(ii, 3), ...
                                      'nearest');
        end
    end
    output_tess_quat = [output_tess_quat_1 output_tess_quat_2 output_tess_quat_3 output_tess_quat_4];
    %
    if(add_local_perturbations == 0)
        local_perturbations = zeros(size(gridpoints_coord, 1), 10);
    else
        % Interpolate perturbations
    end
    %
    % File to write microstructural data
    f = fopen(ms_file_name, 'w');
    gridpoint_counter = 1;
    for ii = 1:size(gridpoints_coord, 1)
        % Write microstructure data in the column format that fwdmodel
        % python code reads.
        % Strain is really right stretch tensor.
        % I am adding a small noise to orientations and strains. Hence
        % the rand() calls.
        fprintf(f, '%12.6f, %12.6f, %12.6f, %s, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f\n', ...
            gridpoints_coord(gridpoint_counter, 1), ...
            gridpoints_coord(gridpoint_counter, 1), ...
            gridpoints_coord(gridpoint_counter, 1), ...
            material_name, ...
            output_tess_quat(gridpoint_counter, 1) + local_perturbations(gridpoint_counter, 1)*sigma_o, ...
            output_tess_quat(gridpoint_counter, 2) + local_perturbations(gridpoint_counter, 2)*sigma_o, ...
            output_tess_quat(gridpoint_counter, 3) + local_perturbations(gridpoint_counter, 3)*sigma_o, ...
            output_tess_quat(gridpoint_counter, 4) + local_perturbations(gridpoint_counter, 4)*sigma_o, ...
            1.0 + strain(output_tess(gridpoint_counter), 1) + local_perturbations(gridpoint_counter, 5)*sigma_e, ...
            1.0 + strain(output_tess(gridpoint_counter), 2) + local_perturbations(gridpoint_counter, 6)*sigma_e, ...
            1.0 + strain(output_tess(gridpoint_counter), 3) + local_perturbations(gridpoint_counter, 7)*sigma_e, ...
            strain(output_tess(gridpoint_counter), 4) + local_perturbations(gridpoint_counter, 8)*sigma_e, ...
            strain(output_tess(gridpoint_counter), 5) + local_perturbations(gridpoint_counter, 9)*sigma_e, ...
            strain(output_tess(gridpoint_counter), 6) + local_perturbations(gridpoint_counter, 10)*sigma_e, ...
            output_tess_intensity_multiplier(gridpoint_counter));
        gridpoint_counter = gridpoint_counter + 1;
    end
    % Don't forget to close files.
    fclose(f);
end
%

%%
% Plot COM, orientations etc? Useful for debugging.
plot_figures = 1;
%
if(plot_figures)
    if(read_grid_data == 0)
        % IPF with orientations
        q = quaternion(quat');
        o = orientation(q, crystalSymmetry('cubic'), specimenSymmetry('1'));
        plotIPDF(o, yvector, 'xAxisDirection','east', 'MarkerSize', 8);
        export_fig 'NiTi-A_orientations' -png -r100
        % Scatter plot with tessellation
        figure;
        scatter3(tess_coords(:, 1), tess_coords(:, 2), tess_coords(:, 3), 4, tess(:), 'filled')
        xlabel('X (um)'); ylabel('y (um)'); zlabel('Z (um)');
        set(gca, 'FontSize', 24);
        axis equal; axis vis3d;
        export_fig 'NiTi-A_tess' -png -r100
        % Scatter plot of the intensity factor
        figure;
        scatter3(tess_coords(:, 1), tess_coords(:, 2), tess_coords(:, 3), 4, intensity_multiplier(:), 'filled')
        colorbar;
        xlabel('X (um)'); ylabel('y (um)'); zlabel('Z (um)');
        set(gca, 'FontSize', 24);
        axis equal; axis vis3d;
        export_fig 'NiTi-A_tess_intensity_mult' -png -r100
        if(add_local_perturbations ~= 0)
            % Scatter plot of unsmoothed perturbations
            figure;
            color_tmp = squeeze(local_perturbations_tmp(:, :, :, 1));
            scatter3(tess_coords(:, 1), tess_coords(:, 2), tess_coords(:, 3), 4, color_tmp(:), 'filled')
            colorbar;
            xlabel('X (um)'); ylabel('y (um)'); zlabel('Z (um)');
            set(gca, 'FontSize', 24);
            axis equal; axis vis3d;
            export_fig 'NiTi-A_tess_perturb' -png -r100
            % Scatter plot of unsmoothed perturbations
            figure;
            color_tmp = squeeze(local_perturbations(:, :, :, 1));
            scatter3(tess_coords(:, 1), tess_coords(:, 2), tess_coords(:, 3), 4, color_tmp(:), 'filled')
            xlabel('X (um)'); ylabel('y (um)'); zlabel('Z (um)');
            set(gca, 'FontSize', 24);
            axis equal; axis vis3d;
            export_fig 'NiTi-A_tess_perturb_smooth' -png -r100
        end
        if(add_powder_to_gb)
            % Scatter plot of gb areas
            figure;
            scatter3(tess_coords(:, 1), tess_coords(:, 2), tess_coords(:, 3), 4, is_gb(:), 'filled')
            xlabel('X (um)'); ylabel('y (um)'); zlabel('Z (um)');
            set(gca, 'FontSize', 24);
            axis equal; axis vis3d;
            export_fig 'NiTi-A_tess_grain_boundary' -png -r100
        end
    else
        % Grid data (which may be irregular) was provided by user.
        % Alternate plotting routine.
        % IPF with orientations
        q = quaternion(quat');
        o = orientation(q, crystalSymmetry('cubic'), specimenSymmetry('1'));
        plotIPDF(o, yvector, 'xAxisDirection','east', 'MarkerSize', 8);
        export_fig 'NiTi-A_orientations' -png -r100
        % Scatter plot with tessellation
        figure;
        scatter3(gridpoints_coord(:, 1), gridpoints_coord(:, 2), gridpoints_coord(:, 3), 4, output_tess, 'filled')
        xlabel('X (um)'); ylabel('y (um)'); zlabel('Z (um)');
        set(gca, 'FontSize', 24);
        axis equal; axis vis3d;
        export_fig 'NiTi-A_tess' -png -r100
        % Scatter plot of the intensity factor
        figure;
        scatter3(gridpoints_coord(:, 1), gridpoints_coord(:, 2), gridpoints_coord(:, 3), 4, output_tess_intensity_multiplier, 'filled')
        colorbar;
        xlabel('X (um)'); ylabel('y (um)'); zlabel('Z (um)');
        set(gca, 'FontSize', 24);
        axis equal; axis vis3d;
        export_fig 'NiTi-A_tess_intensity_mult' -png -r100
        if(add_local_perturbations ~= 0)
            % Scatter plot of unsmoothed perturbations
            figure;
            scatter3(gridpoints_coord(:, 1), gridpoints_coord(:, 2), gridpoints_coord(:, 3), 4, local_perturbations(:, 1), 'filled')
            colorbar;
            xlabel('X (um)'); ylabel('y (um)'); zlabel('Z (um)');
            set(gca, 'FontSize', 24);
            axis equal; axis vis3d;
            export_fig 'NiTi-A_tess_perturb' -png -r100
        end
        if(add_powder_to_gb)
            % Scatter plot of gb areas
            figure;
            scatter3(gridpoints_coord(:, 1), gridpoints_coord(:, 2), gridpoints_coord(:, 3), 4, output_is_gb, 'filled')
            xlabel('X (um)'); ylabel('y (um)'); zlabel('Z (um)');
            set(gca, 'FontSize', 24);
            axis equal; axis vis3d;
            export_fig 'NiTi-A_tess_grain_boundary' -png -r100
        end
    end
end
% Prepare variables for function output
if(~isempty(gridpoints_coord))
    tess_coords = gridpoints_coord;
    tess = output_tess;
    tess_quat = output_tess_quat;
    intensity_multiplier = output_tess_intensity_multiplier;
else
    tess = tess(:);
    tess_quat = reshape(tess_quat, N*N*N, 4);
    strain = reshape(strain, N*N*N, 6);
    intensity_multiplier = intensity_multiplier(:);
end
end
