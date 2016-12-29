% Load/analyse data from MIDAS grain tracking results
%    Read Grains.csv and GrainTracking.csv files for multiple layers and 
%    plot COM, orientations etc. Create a grain correspondence table.
%

parent_grains_file_stem = ['midas_tracked_grains_files/grains_00_layer'];
total_num_layers = 5;
loading_steps = 1:8;
%
graindata_parent = {};
num_grains_all_layers = 0;
num_grains_each_layer = [];
% Read grain data for the parent loading step
for layer_num = 1:total_num_layers
    try
        gfile_name = [parent_grains_file_stem num2str(layer_num) '.csv'];
        grains_layer = load(gfile_name);
        graindata_parent{end + 1} = grains_layer;
        num_grains_all_layers = num_grains_all_layers + size(grains_layer, 1);
        num_grains_each_layer(end + 1) = size(grains_layer, 1);
    catch
        warning('Skipping layer %d', layer_num);
    end
end
num_grains_upto_layer = cumsum(num_grains_each_layer); % Number of grains upto a layer.
%num_grains_upto_layer = num_grains_upto_layer - num_grains_upto_layer(1);
num_grains_upto_layer = [0 num_grains_upto_layer];
num_grains_upto_layer(end) = [];
%
% Create a grain correspondence table
grain_correspondence_table_tracked = zeros(num_grains_all_layers, numel(loading_steps(:)));
grain_correspondence_table_tracked(:, 1) = [1:num_grains_all_layers]';
% For each layer, loop over all loading steps and load grain data and grain
% correspondence data.
num_grains_upto_layer_tracked = zeros(total_num_layers, numel(loading_steps));
grain_tracking_data = cell(total_num_layers, numel(loading_steps));
grain_data = cell(total_num_layers, numel(loading_steps));
%
for layer_num = 1:total_num_layers
    parent_layer_data = graindata_parent{layer_num};
    for ll = 1:numel(loading_steps)
        ll_step = loading_steps(ll);
        grains_file_stem = ['midas_tracked_grains_files/grains_0' num2str(ll_step) '_layer'];
        grains_correspondence_file_stem = ['midas_tracked_grains_files/tracked_grains_0' num2str(ll_step) '_layer'];
        try
            grains_current = load([grains_file_stem num2str(layer_num) '.csv']);
            grain_data{layer_num, ll_step} = grains_current;
            grains_correspondence_current = importdata([grains_correspondence_file_stem num2str(layer_num) '.csv']);
            grains_correspondence_current = grains_correspondence_current.data;
            grain_tracking_data{layer_num, ll_step} = grains_correspondence_current;
            num_grains_upto_layer_tracked(layer_num, ll_step) = size(grains_current, 1);
        catch
            warning('Skipping load step %d layer %d', ll_step, layer_num);
        end
    end
end
%
num_grains_upto_layer_tracked = cumsum(num_grains_upto_layer_tracked, 1);
%for mm = 1:size(num_grains_upto_layer_tracked, 2)
%    num_grains_upto_layer_tracked(:, mm) = num_grains_upto_layer_tracked(:, mm) - num_grains_upto_layer_tracked(1, mm);
%end
for mm = 1:size(num_grains_upto_layer_tracked, 2)
   num_grains_upto_layer_tracked(end, mm) = 0;
end
num_grains_upto_layer_tracked = circshift(num_grains_upto_layer_tracked, 1, 1);
% Now we have all the data needed to create a grain correspondence table
for layer_num = 1:total_num_layers
    for ll = 1:numel(loading_steps)
        % Get parent and tracked grain data
        parent_data_ll = graindata_parent{layer_num};
        current_data_ll = grain_data{layer_num, ll};
        % Get the tracking data (old_grain_id <-> new_grain_id)
        tracking_data_ll = grain_tracking_data{layer_num, ll};
        % Loop over all correspondences
        for mm = 1:size(tracking_data_ll, 1)
            % Non-zero new grain id means that the grain was tracked
            if(tracking_data_ll(mm, 2) ~= 0)
                parent_grain_id = tracking_data_ll(mm, 1);
                current_grain_id = tracking_data_ll(mm, 2);
                % Matching grain index in the layer-specific numbering
                p = find(parent_data_ll(:, 1) == parent_grain_id);
                t = find(current_data_ll(:, 1) == current_grain_id);
                % Transform the ids to the global numbering for that load
                p = p + num_grains_upto_layer(layer_num);
                t = t + num_grains_upto_layer_tracked(layer_num, ll);
                grain_correspondence_table_tracked(p, ll + 1) = t;
            end
        end
    end
end
%
generate_tracked_grain_data = 1;
if(generate_tracked_grain_data)
    APS_2016_s1_cycle_1_radius_grains = {};
    APS_2016_s1_cycle_1_lattice_strain = [];
    APS_2016_s1_cycle_1_lattice_strain_surf = [];
    APS_2016_s1_cycle_1_lattice_strain_interior = [];
    APS_2016_s1_cycle_1_lattice_strain_stdev = [];
    APS_2016_s1_cycle_1_lattice_strain_stdev_surf = [];
    APS_2016_s1_cycle_1_lattice_strain_stdev_interior = [];
    APS_2016_s1_cycle_1_num_of_grains = [];
    APS_2016_s1_cycle_1_num_of_grains_surf = [];
    APS_2016_s1_cycle_1_orientations = {};
    APS_2016_s1_cycle_1_com = {};
    APS_2016_s1_cycle_1_lattice_strain_grains = {};
    APS_2016_s1_cycle_1_lattice_strain_full_grains = {};
    APS_2016_s1_cycle_1_grain_id_correspondence = {};

    for ll = [0:8]
        % INPUTS
        % Load stiffness tensor for NiTi from file
        stiffness_tensor = loadTensor('NiTi_cubic_elastic_constants.data' , crystalSymmetry('cubic'), 'name', 'NiTi stiffness');
        total_num_layers = 5;
        layer_thickness  = 150;
        grains_file_stem = ['midas_tracked_grains_files/grains_0' num2str(ll) '_layer'];                 % grains file name = stem + layer_num
        grains_to_plot = [];                                                   % Array to plot specific grains, empty to plot all
        make_plots = false;
        % Read Grains.csv data
        grains = [];
        grains_layer_id = [];
        grains_strains_fab_matrix = [];
        for layer_num = 1:total_num_layers
            try
                grains_layer = load([grains_file_stem num2str(layer_num) '.csv']);
                grains_layer(:, 13) = grains_layer(:, 13) + layer_thickness*(layer_num - 1);
                grains = [grains; grains_layer];
                %grain_id{end + 1, layer_num} = grains(:, 1);
                grains_layer_id = [grains_layer_id; layer_num * ones(size(grains_layer, 1), 1)];
            catch
                warning('Skipping load step %d layer %d', ll, layer_num);
            end
        end
        % Store useful quantities in arrays
        % Note: after applying the RESRF2APS transformation, loading axis is Y,
        % just like heXRD
        RESRF2APS   = RMatOfQuat(QuatOfESRF2APS);
        %
        num_grains     = size(grains, 1);
        %
        grains_com = RESRF2APS*grains(:, 11:13)';
        grains_com = grains_com';
        grains_lattice_params = grains(:, 14:19);
        grains_strains_fab = grains(:, 25:33);
        % After applying ESRF (FABLE) -> APS transformation, grain strain
        % tensor in the laboratory coordinates. Y (2, 2) is loading axis.
        for i = 1:1:num_grains
            grains_strains_fab_matrix(:,:,i)   =  RESRF2APS*reshape(grains_strains_fab(i, :), 3, 3)'*RESRF2APS';
        end
        grains_strains_ken = grains(:, 34:42);
        grains_orient = grains(:, 2:10);
        grains_radius = grains(:, 23);
        grains_completeness = grains(:, 24);
        grains_radial_dist = sqrt(grains_com(:, 1).^2 + grains_com(:, 3).^2);

        disp(' ')
        disp(['Loading step ' num2str(ll)])
        disp(['Total number of grains = ' num2str(size(grains, 1))])
        disp(['Grains with radial dist > 300 = ' num2str(sum(grains_radial_dist > 300))])
        disp(['Fraction of surface grains = ' num2str(sum(grains_radial_dist > 300)/size(grains, 1))])
        disp(['Mean radius = ' num2str(mean(grains_radius, 1))])
        disp(['Mean radius (surf) = ' num2str(mean(grains_radius(grains_radial_dist(:) > 300), 1))])
        disp(['Mean strains (FABLE XYZ, Y = loading) = ' num2str(mean(squeeze(grains_strains_fab_matrix(1, 1, :)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(2, 2, :)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(3, 3, :))))])
        disp(['Mean strains (surface) (FABLE XYZ, Y = loading) = ' num2str(mean(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300))))])
        disp(['Mean strains (interior) (FABLE XYZ, Y = loading) = ' num2str(mean(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) < 300)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) < 300)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) < 300))))])
        disp(['standard dev strains (FABLE XYZ, Y = loading) = ' num2str(std(squeeze(grains_strains_fab_matrix(1, 1, :)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(2, 2, :)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(3, 3, :))))])
        disp(['standard dev strains (surface) (FABLE XYZ, Y = loading) = ' num2str(std(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300))))])
        disp(['standard dev strains (interior) (FABLE XYZ, Y = loading) = ' num2str(std(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) < 300)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) < 300)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) < 300))))])
        %disp(['Mean strains (Kenesei) = ' num2str(mean(grains_strains_ken(:, 1), 1)) ', ' num2str(mean(grains_strains_ken(:, 2), 1)) ', ' num2str(mean(grains_strains_ken(:, 3), 1))])
        disp(' ')
        % Save data
        APS_2016_s1_cycle_1_lattice_strain(end + 1, :) = [mean(squeeze(grains_strains_fab_matrix(1, 1, :))) mean(squeeze(grains_strains_fab_matrix(2, 2, :))) mean(squeeze(grains_strains_fab_matrix(3, 3, :)))];
        APS_2016_s1_cycle_1_lattice_strain_surf(end +1, :) = [mean(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300))) mean(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300))) mean(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300)))];
        APS_2016_s1_cycle_1_lattice_strain_interior(end +1, :) = [mean(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) < 300))) mean(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) < 300))) mean(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) < 300)))];
        APS_2016_s1_cycle_1_lattice_strain_grains{end + 1} = [squeeze(grains_strains_fab_matrix(1, 1, :)) squeeze(grains_strains_fab_matrix(2, 2, :)) squeeze(grains_strains_fab_matrix(3, 3, :))];
        APS_2016_s1_cycle_1_lattice_strain_full_grains{end + 1} = grains_strains_fab_matrix;
        APS_2016_s1_cycle_1_lattice_strain_stdev(end + 1, :) = [std(squeeze(grains_strains_fab_matrix(1, 1, :))) std(squeeze(grains_strains_fab_matrix(2, 2, :))) std(squeeze(grains_strains_fab_matrix(3, 3, :)))];
        APS_2016_s1_cycle_1_lattice_strain_stdev_surf(end +1, :) = [std(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300))) std(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300))) std(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300)))];
        APS_2016_s1_cycle_1_lattice_strain_stdev_interior(end +1, :) = [std(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) < 300))) std(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) < 300))) std(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) < 300)))];
        APS_2016_s1_cycle_1_num_of_grains(end +1) = size(grains, 1);
        APS_2016_s1_cycle_1_num_of_grains_surf(end + 1) = sum(grains_radial_dist > 300);
        APS_2016_s1_cycle_1_com{end + 1} = grains_com;
        APS_2016_s1_cycle_1_radius_grains{end + 1} = grains_radius;

        % Calculate the orientation (MTEX) for each grain
        RMats = [];
        for i = 1:1:num_grains
            RMats(:,:,i)   =  RESRF2APS*reshape(grains(i,2:10), 3, 3)';
        end
    %     qsym    = CubSymmetries; Rsym    = RMatOfQuat(qsym);
    %     quat    = ToFundamentalRegionQ(QuatOfRMat(RMats), qsym);
        quat    = QuatOfRMat(RMats);
        cs = crystalSymmetry('m-3m');
        ss = specimenSymmetry('triclinic');
        q = quaternion(quat);
        r = rotation(q);
        o = orientation(q, cs, ss);
        APS_2016_s1_cycle_1_orientations{end + 1} = o;
        %
        oM = ipdfHSVOrientationMapping(crystalSymmetry('cubic'));
        % Colors are generated for [0 0 1] IPF. Need to fix that.
        oM.inversePoleFigureDirection = vector3d(0,1,0);
        rgb = oM.orientation2color(o);
        % Rotate stiffness tensor to global coords in each gran
        % Calculate stiffness along loading direction (in MPa)
    %     for i = 1:num_grains
    %         C_rot = rotate(stiffness_tensor, r(i));
    %         stiffness_loading(i) = YoungsModulus(C_rot , yvector)/1e6;
    %     end
        if(make_plots)
            % Visualize IPF
            %plotIPDF(o, yvector, 'Property', stiffness_loading'.*squeeze(grains_strains_fab_matrix(2, 2, :)), 'xAxisDirection','east', 'MarkerSize', 4);
            plotIPDF(o, yvector, 'Property', squeeze(grains_strains_fab_matrix(2, 2, :)), 'xAxisDirection','east', 'MarkerSize', 4);
            colorbar; colormap jet; caxis([-3000 5000]);
            print(gcf, ['ipf_strain_step_' num2str(ll)], '-dpng');
            % Visualize COM
            figure;
            scatter3(grains_com(:, 1), grains_com(:, 2), grains_com(:, 3), 121, grains_radial_dist > 400, 'filled');
            campos([0 1800 0]); camtarget([0 200 0]); camup([0 0 1]);
            axis equal;
            axis vis3d;
        end
        close all;
        %clear q r o RMats grains_strains_fab_matrix;
    end
end

