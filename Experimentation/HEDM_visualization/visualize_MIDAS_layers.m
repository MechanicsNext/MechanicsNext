% Visualize MIDAS grain reconstruction for multiple layers + loading steps
%
%    This script assumes that ff-HEDM measurement was performed for
%    multiple non-overlapping layers at each loading step and the
%    measurements were repeated at multiple loading steps.
%    ________             ________             ________
%   |________| Layer 4   |________| Layer 4   |________| Layer 4
%   |________| Layer 3   |________| Layer 3   |________| Layer 3
%   |________| Layer 2   |________| Layer 2   |________| Layer 2
%   |________| Layer 1   |________| Layer 1   |________| Layer 1
%     Loading step 0  ->   Loading step 1  ->   Loading step 2  -> ...
%
% PREREQUISITES
%   MTEX (https://mtex-toolbox.github.io/)
%
% INPUTS
grains_file_format      = 'grains_files/Grains_Step_%d_Layer%d.csv';       % Format of the grains file name (including the absolute path). See details below.
number_of_loading_steps = 9;                                               % Total number of loading steps
first_loading_step      = 0;                                               % Loading steps will go from first to first + number_of_loading_steps - 1
num_layers_each_step    = 5;                                               % Number of layers for each load step
first_layer             = 1;                                               % Number of the first layer at a given loading step.
layer_thickness         = 150;                                             % Thickness of a lyer in microns
space_group             = 'cubic';                                         % Space group of the material in MTEX format. e.g., 'cubic', '11', 'hexagonal'
%
% 'grains_file_format' must be in a format that can be used with sprintf
% command. sprintf will be called with two integers
% sprintf(grains_file_format, current_load_num, current_layer_num). Last
% two arguments are integers corresponding to the loading step and layer
% number.
% For example, if the files are named Grains_Step_1_layer_1.csv,
% Grains_Step_1_layer_2.csv etc. with loading steps 0 to 5 and layers 
% 1 to 3 then:
% grains_file_format = 'Grains_Step_%d_layer_%d.csv';
% number_of_loading_steps = 6;
% first_loading_step = 0;
% num_layers_each_step = 3;
% first_layer = 1;
%
% Initialize storage arrays
aggregated_grain_radius_grains = {};
aggregated_grain_lattice_strain = [];
aggregated_grain_lattice_strain_surf = [];
aggregated_grain_lattice_strain_interior = [];
aggregated_grain_lattice_strain_stdev = [];
aggregated_grain_lattice_strain_stdev_surf = [];
aggregated_grain_lattice_strain_stdev_interior = [];
aggregated_grain_num_of_grains = [];
aggregated_grain_num_of_grains_surf = [];
aggregated_grain_orientations = {};
aggregated_grain_com = {};
aggregated_grain_lattice_strain_grains = {};
aggregated_grain_lattice_strain_full_grains = {};

for ll = [first_loading_step:(number_of_loading_steps + first_loading_step - 1)]
    % INPUTS
    total_num_layers = num_layers_each_step;
    grains_to_plot = [];                                                       % Array to plot specific grains, empty to plot all
    make_plots = false;
    % Read Grains.csv data
    grains = [];
    grains_layer_id = [];
    grains_strains_fab_matrix = [];
    for layer_num = first_layer:(total_num_layers + first_layer - 1)
        grains_layer = load(sprintf(grains_file_format, ll, layer_num));
        grains_layer(:, 13) = grains_layer(:, 13) + layer_thickness*(layer_num - 1);
        grains = [grains; grains_layer];
        grains_layer_id = [grains_layer_id; layer_num * ones(size(grains_layer, 1), 1)];
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
    aggregated_grain_lattice_strain(end + 1, :) = [mean(squeeze(grains_strains_fab_matrix(1, 1, :))) mean(squeeze(grains_strains_fab_matrix(2, 2, :))) mean(squeeze(grains_strains_fab_matrix(3, 3, :)))];
    aggregated_grain_lattice_strain_surf(end +1, :) = [mean(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300))) mean(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300))) mean(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300)))];
    aggregated_grain_lattice_strain_interior(end +1, :) = [mean(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) < 300))) mean(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) < 300))) mean(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) < 300)))];
    aggregated_grain_lattice_strain_grains{end + 1} = [squeeze(grains_strains_fab_matrix(1, 1, :)) squeeze(grains_strains_fab_matrix(2, 2, :)) squeeze(grains_strains_fab_matrix(3, 3, :))];
    aggregated_grain_lattice_strain_full_grains{end + 1} = grains_strains_fab_matrix;
    aggregated_grain_lattice_strain_stdev(end + 1, :) = [std(squeeze(grains_strains_fab_matrix(1, 1, :))) std(squeeze(grains_strains_fab_matrix(2, 2, :))) std(squeeze(grains_strains_fab_matrix(3, 3, :)))];
    aggregated_grain_lattice_strain_stdev_surf(end +1, :) = [std(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300))) std(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300))) std(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300)))];
    aggregated_grain_lattice_strain_stdev_interior(end +1, :) = [std(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) < 300))) std(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) < 300))) std(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) < 300)))];
    aggregated_grain_num_of_grains(end +1) = size(grains, 1);
    aggregated_grain_num_of_grains_surf(end + 1) = sum(grains_radial_dist > 300);
    aggregated_grain_com{end + 1} = grains_com;
    aggregated_grain_radius_grains{end + 1} = grains_radius;

    % Calculate the orientation (MTEX) for each grain
    RMats = [];
    for i = 1:1:num_grains
        RMats(:,:,i)   =  RESRF2APS*reshape(grains(i,2:10), 3, 3)';
    end
    %
    quat = QuatOfRMat(RMats);
    cs = crystalSymmetry(space_group);
    ss = specimenSymmetry('1');
    q = quaternion(quat);
    r = rotation(q);
    o = orientation(q, cs, ss);
    aggregated_grain_orientations{end + 1} = o;
    %
    oM = ipdfHSVOrientationMapping(crystalSymmetry('432'));
    % Colors are generated for [0 0 1] IPF. Need to fix that.
    oM.inversePoleFigureDirection = vector3d(0,1,0);
    rgb = oM.orientation2color(o);
    if(make_plots)
        % Visualize IPF with axial (Y) lattice strain as color
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
end

